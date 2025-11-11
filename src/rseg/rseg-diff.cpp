/*
 * Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song and Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

static constexpr auto about = R"(rseg-diff

Segment the genome according to difference in mapped read density.
)";

#include "EvaluateBoundaries.hpp"
#include "GenomicRegion.hpp"
#include "LoadReadsByRegion.hpp"
#include "OptionParser.hpp"

#include "ReadCounts.hpp"
#include "SelectBinSize.hpp"
#include "SplitDistro.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include "ModelParams.hpp"

#include "TwoStateScaleSplitHMM.hpp"
#include "TwoStateScaleSplitResolveMixture.hpp"

#include "ThreeStateScaleSplitHMM.hpp"
#include "ThreeStateScaleSplitResolveMixture.hpp"

#include "rseg_utils.hpp"

#include "CLI11.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <print>
#include <random>
#include <utility>

// functions for two-state modes

void
output_boundaries(
  const std::vector<double> &tmp_read_bins, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points,
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const TwoStateScaleSplitHMM &hmm, const std::vector<SplitDistro> &distros,
  const std::vector<double> &start_trans,
  const std::vector<std::vector<double>> &trans,
  const std::vector<double> &end_trans, const std::vector<bool> &tmp_classes,
  const std::string &boundary_file, const std::string &boundary_score_file,
  const bool VERBOSE) {

  std::vector<std::vector<std::vector<double>>> post_trans_scores;
  hmm.TransitionPosteriors(tmp_read_bins, scales, reset_points, start_trans,
                           trans, end_trans, distros.front(), distros.back(),
                           post_trans_scores);

  std::vector<double> &f_to_f_scores(post_trans_scores[0][0]);
  std::vector<double> &f_to_b_scores(post_trans_scores[0][1]);
  std::vector<double> &b_to_f_scores(post_trans_scores[1][0]);
  std::vector<double> &b_to_b_scores(post_trans_scores[1][1]);

  std::vector<double> tmp_boundary_scores(f_to_b_scores.size());
  std::vector<int> transitions;
  for (size_t i = 0; i < tmp_classes.size(); ++i)
    if (i == 0 || tmp_classes[i] != tmp_classes[i - 1])
      transitions.push_back(i);

  transitions.push_back(tmp_classes.size());
  size_t j = 0;

  for (int i = 0; static_cast<size_t>(i) < f_to_b_scores.size(); ++i) {
    static const size_t sz = tmp_classes.size();

    if (abs(i - transitions[j]) > abs(i - transitions[j + 1]))
      ++j;
    const size_t trn = transitions[j];
    if (trn == 0 || trn == sz)
      tmp_boundary_scores[i] = 0;
    else
      tmp_boundary_scores[i] =
        (tmp_classes[trn]) ? b_to_f_scores[i] : f_to_b_scores[i];
  }

  std::vector<std::vector<double>> boundary_scores;
  expand_bins(tmp_boundary_scores, reset_points, boundary_scores);

  //// generate control sample
  const size_t rand_sample_size =
    std::min(size_t(150000), tmp_read_bins.size());
  std::vector<std::pair<double, double>> vals_control(rand_sample_size);

  for (size_t i = 0; i < rand_sample_size; ++i)
    vals_control[i] = std::make_pair(tmp_read_bins[i], scales[i]);

  std::vector<size_t> reset_points_control;
  for (size_t i = 0; i < reset_points.size(); ++i)
    if (reset_points[i] < rand_sample_size)
      reset_points_control.push_back(reset_points[i]);
    else
      break;
  reset_points_control.push_back(rand_sample_size);

  std::random_device rd{};
  std::mt19937 g{rd()};

  for (size_t i = 0; i < reset_points_control.size() - 1; ++i)
    std::shuffle(vals_control.begin() + reset_points_control[i],
                 vals_control.begin() + reset_points_control[i + 1], g);

  std::vector<double> read_counts_control(rand_sample_size),
    scales_control(rand_sample_size);

  for (size_t i = 0; i < rand_sample_size; ++i) {
    read_counts_control[i] = vals_control[i].first;
    scales_control[i] = vals_control[i].second;
  }
  vals_control.clear();

  std::vector<std::vector<std::vector<double>>> post_trans_scores_control;
  hmm.TransitionPosteriors(read_counts_control, scales_control,
                           reset_points_control, start_trans, trans, end_trans,
                           distros.front(), distros.back(),
                           post_trans_scores_control);

  std::vector<double> &f_to_b_scores_control = post_trans_scores_control[0][1];
  std::vector<double> &b_to_f_scores_control = post_trans_scores_control[1][0];

  std::vector<double> boundary_scores_control(f_to_b_scores_control.size(), 0);

  for (size_t i = 0; i < boundary_scores_control.size(); ++i)
    boundary_scores_control[i] =
      std::max(f_to_b_scores_control[i], b_to_f_scores_control[i]);

  std::sort(boundary_scores_control.begin(), boundary_scores_control.end());
  double fdr = 0.05;

  double cutoff = boundary_scores_control[static_cast<size_t>(
    boundary_scores_control.size() * (1 - fdr))];

  read_counts_control.clear();
  scales_control.clear();
  b_to_f_scores_control.clear();
  f_to_b_scores_control.clear();
  boundary_scores_control.clear();
  post_trans_scores_control.clear();

  //// finish generating control sample

  std::vector<GenomicRegion> boundaries;
  BoundEval be(1, 1);
  be.evaluate(bin_bounds, reset_points, tmp_classes, tmp_boundary_scores,
              f_to_f_scores, f_to_b_scores, b_to_f_scores, b_to_b_scores,
              cutoff, true, boundaries);

  // write result files
  WriteBEDFile(boundary_file, boundaries);
  if (VERBOSE)
    std::cout << "Boundary file: " + boundary_file << '\n';

  if (!boundary_score_file.empty() && boundary_score_file != "None") {
    write_wigfile(boundary_scores, bin_bounds, boundary_score_file);
    if (VERBOSE)
      std::cout << "Boundary score file: " + boundary_score_file << '\n';
  }
}

void
output_domains(const std::vector<double> &tmp_read_bins,
               const std::vector<double> &tmp_scales,
               const std::vector<size_t> &reset_points,
               const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
               const TwoStateScaleSplitHMM &hmm,
               const std::vector<SplitDistro> &distros,
               const std::vector<double> &start_trans,
               const std::vector<std::vector<double>> &trans,
               const std::vector<double> &end_trans,
               const std::vector<bool> &tmp_classes,
               const double posterior_cutoff, const size_t undef_region_cutoff,
               const double cdf_cutoff, const std::string &domain_file,
               const std::string &posterior_score_file, const bool VERBOSE) {
  // Obtain the scores for the current domain class
  std::vector<double> tmp_scores;
  hmm.PosteriorScores(tmp_read_bins, tmp_scales, reset_points, start_trans,
                      trans, end_trans, distros.front(), distros.back(),
                      tmp_classes, tmp_scores);

  std::vector<std::vector<double>> read_bins, scores, scales;
  std::vector<std::vector<bool>> classes;
  expand_bins(tmp_read_bins, reset_points, read_bins);
  expand_bins(tmp_scores, reset_points, scores);
  expand_bins(tmp_scales, reset_points, scales);
  expand_bins(tmp_classes, reset_points, classes);

  std::vector<std::vector<GenomicRegion>> domains;
  build_domains(bin_bounds, classes, scores, posterior_cutoff, domains,
                undef_region_cutoff);
  pick_domains(bin_bounds, read_bins, scales, distros, domains, cdf_cutoff);

  // output domains
  write_bed_file(domains, domain_file);
  if (VERBOSE)
    std::cout << "Domains file: " + domain_file << '\n';

  if (!posterior_score_file.empty() && posterior_score_file != "None") {
    size_t k = 0;
    for (size_t i = 0; i < scores.size(); ++i)
      for (size_t j = 0; j < scores[i].size(); ++j) {
        scores[i][j] == tmp_classes[k] ? scores[i][j] : 1 - scores[i][j];
        ++k;
      }
    write_wigfile(scores, bin_bounds, posterior_score_file);
    if (VERBOSE)
      std::cout << "Bin score file: " + posterior_score_file << '\n';
  }
}

// end of functions for two-state modes

// for three-state mode
void
output_boundaries(
  const std::vector<double> &tmp_read_bins, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points,
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const ThreeStateScaleSplitHMM &hmm, const std::vector<SplitDistro> &distros,
  const std::vector<double> &start_trans,
  const std::vector<std::vector<double>> &trans,
  const std::vector<double> &end_trans, const std::vector<size_t> &classes,
  const std::string &boundary_file, const std::string &boundary_score_file,
  const bool VERBOSE) {

  const size_t NUM_OF_STATES = 3;

  std::vector<std::vector<std::vector<double>>> post_trans_scores;
  hmm.TransitionPosteriors(tmp_read_bins, scales, reset_points, start_trans,
                           trans, end_trans, distros.front(), distros[1],
                           distros.back(), post_trans_scores);

  std::vector<size_t> change_points;
  for (size_t i = 0; i < classes.size(); ++i)
    if (i == 0 || classes[i] != classes[i - 1])
      change_points.push_back(i);
  change_points.push_back(classes.size());

  std::vector<double> tmp_boundary_scores(tmp_read_bins.size());
  size_t j = 0;

  for (int i = 0; static_cast<size_t>(i) < tmp_boundary_scores.size(); ++i) {
    static const size_t sz = tmp_boundary_scores.size();

    if (abs(i - change_points[j]) > abs(i - change_points[j + 1]))
      ++j;

    const size_t trn = change_points[j];

    if (trn == 0 || trn == sz)  // indicates a starting domain
      tmp_boundary_scores[i] = 0;
    else
      tmp_boundary_scores[i] =
        post_trans_scores[classes[trn - 1]][classes[trn]][i];
  }

  std::vector<std::vector<double>> boundary_scores;
  expand_bins(tmp_boundary_scores, reset_points, boundary_scores);

  ///// genereate contronl sample
  const size_t rand_sample_size =
    std::min(size_t(150000), tmp_read_bins.size());
  std::vector<std::pair<double, double>> vals_control(rand_sample_size);

  for (size_t i = 0; i < rand_sample_size; ++i)
    vals_control[i] = std::make_pair(tmp_read_bins[i], scales[i]);

  std::vector<size_t> reset_points_control;
  for (size_t i = 0; i < reset_points.size(); ++i)
    if (reset_points[i] < rand_sample_size)
      reset_points_control.push_back(reset_points[i]);
    else
      break;
  reset_points_control.push_back(rand_sample_size);

  std::random_device rd{};
  std::mt19937 g{rd()};

  for (size_t i = 0; i < reset_points_control.size() - 1; ++i)
    std::shuffle(vals_control.begin() + reset_points_control[i],
                 vals_control.begin() + reset_points_control[i + 1], g);

  std::vector<double> read_counts_control(rand_sample_size),
    scales_control(rand_sample_size);

  for (size_t i = 0; i < rand_sample_size; ++i) {
    read_counts_control[i] = vals_control[i].first;
    scales_control[i] = vals_control[i].second;
  }
  vals_control.clear();

  std::vector<std::vector<std::vector<double>>> post_trans_scores_control;
  hmm.TransitionPosteriors(read_counts_control, scales_control,
                           reset_points_control, start_trans, trans, end_trans,
                           distros.front(), distros[1], distros.back(),
                           post_trans_scores_control);

  std::vector<double> boundary_scores_control(read_counts_control.size(), 0);
  for (size_t i = 0; i < boundary_scores_control.size(); ++i)
    for (size_t j = 0; j < NUM_OF_STATES; ++j)
      for (size_t k = 0; k < NUM_OF_STATES; ++k)
        if (j != k &&
            post_trans_scores_control[j][k][i] > boundary_scores_control[i])
          boundary_scores_control[i] = post_trans_scores_control[j][k][i];

  std::sort(boundary_scores_control.begin(), boundary_scores_control.end());
  double fdr = 0.05;
  double cutoff = boundary_scores_control[static_cast<size_t>(
    boundary_scores_control.size() * (1 - fdr))];

  read_counts_control.clear();
  scales_control.clear();
  reset_points_control.clear();
  post_trans_scores_control.clear();
  boundary_scores_control.clear();
  //// finish generating control sample

  std::vector<GenomicRegion> boundaries;
  BoundEval be(1, 1);
  be.evaluate(bin_bounds, reset_points, classes, tmp_boundary_scores,
              post_trans_scores, cutoff, boundaries);

  // write result files
  WriteBEDFile(boundary_file, boundaries);
  if (VERBOSE)
    std::cout << "Boundary file: " + boundary_file << '\n';

  if (!boundary_score_file.empty() && boundary_score_file != "None") {
    write_wigfile(boundary_scores, bin_bounds, boundary_score_file);
    if (VERBOSE)
      std::cout << "Boundary score file: " + boundary_score_file << '\n';
  }
}

void
output_domains(const std::vector<double> &tmp_read_bins,
               const std::vector<double> &tmp_scales,
               const std::vector<size_t> &reset_points,
               const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
               const ThreeStateScaleSplitHMM &hmm,
               const std::vector<SplitDistro> &distros,
               const std::vector<double> &start_trans,
               const std::vector<std::vector<double>> &trans,
               const std::vector<double> &end_trans,
               const std::vector<size_t> &classes,
               const double posterior_cutoff, const size_t undef_region_cutoff,
               const double cdf_cutoff, const std::string &domain_file,
               const std::string &posterior_score_file, const bool VERBOSE) {
  // Obtain the scores for the current domain class
  std::vector<double> tmp_scores;
  hmm.PosteriorScores(tmp_read_bins, tmp_scales, reset_points, start_trans,
                      trans, end_trans, distros.front(), distros[1],
                      distros.back(), classes, tmp_scores);

  std::vector<std::vector<double>> read_bins, scores, scales;
  std::vector<std::vector<size_t>> expanded_classes;
  expand_bins(tmp_read_bins, reset_points, read_bins);
  expand_bins(tmp_scores, reset_points, scores);
  expand_bins(tmp_scales, reset_points, scales);
  expand_bins(classes, reset_points, expanded_classes);

  std::vector<std::vector<GenomicRegion>> domains;
  build_domains(bin_bounds, expanded_classes, scores, posterior_cutoff, domains,
                undef_region_cutoff);
  pick_domains_3s(bin_bounds, read_bins, scales, distros, domains, cdf_cutoff);

  // output domains
  write_bed_file(domains, domain_file);
  if (VERBOSE)
    std::cout << "Domains file: " + domain_file << '\n';

  if (!posterior_score_file.empty() && posterior_score_file != "None") {
    std::vector<double> fg_scores, bg_scores;
    hmm.PosteriorScores(tmp_read_bins, tmp_scales, reset_points, start_trans,
                        trans, end_trans, distros.front(), distros[1],
                        distros.back(), fg_scores, bg_scores);
    write_wigfile(fg_scores, bg_scores, bin_bounds, posterior_score_file);
    if (VERBOSE)
      std::cout << "Bin score file: " + posterior_score_file << '\n';
  }
}

// end of functions for three-state modes
int
main(int argc, char *argv[]) {
  static constexpr auto usage = "Usage: rseg-diff [options]";
  // names of emission distributions to use
  static constexpr auto fg_name = "nbdiff";
  static constexpr auto bg_name = "nbdiff";
  static constexpr auto both_domain_ends{true};

  try {

    std::string deads_file;
    std::string chroms_file;
    std::string in_param_file;
    std::string out_param_file;

    // expected size of a domain
    double fg_size = 20000;

    // flags
    bool USE_POSTERIOR{};
    bool REMOVE_JACKPOT{true};
    bool VERBOSE{};
    bool BAM_FORMAT{};

    std::string domain_file;
    std::string posterior_score_file;
    std::string boundary_file;
    std::string boundary_score_file;
    std::string read_counts_file;

    // mode
    int mode = 2;
    const int TEST_CONTROL_MODE = 2;
    const int TEST_TEST_MODE = 3;

    // name of emission distributions
    std::string fg_name{fg_name};
    std::string bg_name{bg_name};

    size_t desert_size = 20000;
    size_t bin_size_step = 50;
    size_t bin_size = 0;
    size_t FRAGMENT_LEN = 0;
    bool waterman = false;
    bool hideaki = false;
    bool hideaki_emp = false;
    bool smooth = true;
    size_t max_iterations = 20;
    size_t training_size = 0;
    double tolerance = 1e-20;
    double min_prob = 1e-20;

    double max_dead_proportion = 0.5;

    double posterior_cutoff = 0.5;
    size_t undef_region_cutoff = 3000;
    double cdf_cutoff = 0.1;

    // Determines how many iterations are used during the initialization
    // phase to find good starting values for the HMM
    const size_t MAX_INITIALIZATION_ITR = 5;

    CLI::App app{about};
    argv = app.ensure_utf8(argv);
    app.usage(usage);

    // clang-format off
    app.set_help_flag("-h,--help", "Print a detailed help message and exit");
    app.add_option("-o,--out", domain_file, "output file")
      ->required();
    app.add_option("--score", posterior_score_file, "Posterior scores file");
    app.add_option("--readcount", read_counts_file, "readcounts file");
    app.add_option("--boundary", boundary_file, "domain boundary file");
    app.add_option("--boundary-score", boundary_score_file, "boundary transition scores file");
    app.add_flag("-v,--verbose", VERBOSE, "print more info");
    // clang-format on

    if (argc == 1) {
      std::println("{}", app.help());
      return EXIT_FAILURE;
    }

    CLI11_PARSE(app, argc, argv);

    ////////////////////// COMMAND LINE OPTIONS /////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "segment the genome according to differential "
                           "mapped read density",
                           "<mapped-read-locations-A> "
                           "<mapped-read-locations-B>");
    opt_parse.add_opt("out", 'o', "domain output file", false, domain_file);
    opt_parse.add_opt("score", '\0', "Posterior scores file", false,
                      posterior_score_file);
    opt_parse.add_opt("readcount", '\0', "readcounts file", false,
                      read_counts_file);
    opt_parse.add_opt("boundary", '\0', "domain boundary file", false,
                      boundary_file);
    opt_parse.add_opt("boundary-score", '\0', "boundary transition scores file",
                      false, boundary_score_file);
    opt_parse.add_opt("chrom", 'c', "file with chromosome sizes (BED format)",
                      true, chroms_file);
    opt_parse.add_opt("deadzones", 'd', "file of deadzones (BED format)", false,
                      deads_file);
    opt_parse.add_opt("bam", 'B', "Input reads file is BAM format", false,
                      BAM_FORMAT);
    opt_parse.add_opt("param-in", '\0', "Input parameters file", false,
                      in_param_file);
    opt_parse.add_opt("param-out", '\0', "Output parameters file", false,
                      out_param_file);
    opt_parse.add_opt("mode", 'm', "running mode 2:test-control; 3: test-test",
                      false, mode);
    opt_parse.add_opt("maxitr", 'i', "maximum iterations for training", false,
                      max_iterations);
    opt_parse.add_opt("bin-size", 'b', "bin size (default: based on data)",
                      false, bin_size);
    opt_parse.add_opt("bin-step", '\0',
                      "minimum bin size (default: " + toa(bin_size_step) + ")",
                      false, bin_size_step);
    opt_parse.add_opt("duplicates", '\0', "keep duplicate reads", false,
                      REMOVE_JACKPOT);
    opt_parse.add_opt("fragment_length", '\0',
                      "Extend reads to fragment length (default not to extend)",
                      false, FRAGMENT_LEN);
    opt_parse.add_opt("Waterman", '\0', "use Waterman's method for bin size",
                      false, waterman);
    opt_parse.add_opt("Hideaki", '\0', "use Hideaki's method for bin size",
                      false, hideaki);
    opt_parse.add_opt("Hideaki-emp", '\0',
                      "use Hideaki's empirical method (default)", false,
                      hideaki_emp);
    opt_parse.add_opt("smooth", '\0',
                      "Indicate whether the rate curve is assumed smooth",
                      false, smooth);
    opt_parse.add_opt("max-dead", '\0',
                      "max deadzone proportion for retained bins", false,
                      max_dead_proportion);
    opt_parse.add_opt("domain-size", 's',
                      "expected domain size "
                      "(default: " +
                        toa(fg_size) + ")",
                      false, fg_size);
    opt_parse.add_opt("desert", 'S',
                      "desert size "
                      "(default: " +
                        toa(desert_size) + ")",
                      false, desert_size);
    opt_parse.add_opt("fg", 'F', "foreground emission distribution", false,
                      fg_name);
    opt_parse.add_opt("bg", 'B', "background emission distribution", false,
                      bg_name);
    opt_parse.add_opt("training-size", '\0',
                      "Max number of data points for training (default: all)",
                      false, training_size);
    opt_parse.add_opt("posterior", 'P',
                      "use posterior decoding "
                      "(default: Viterbi)",
                      false, USE_POSTERIOR);
    opt_parse.add_opt("posterior-cutoff", '\0',
                      "Posterior threshold for signigicant bins", false,
                      posterior_cutoff);
    opt_parse.add_opt("undefined", '\0', "min size of unmappable region", false,
                      undef_region_cutoff);
    opt_parse.add_opt("cutoff", '\0', "cutoff in cdf for identified domains",
                      false, cdf_cutoff);
    opt_parse.add_opt("verbose", 'v', "print more run information", false,
                      VERBOSE);

    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);

    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }

    if (leftover_args.size() < 2) {
      std::cerr << "Need two reads files" << '\n';
      return EXIT_SUCCESS;
    }

    const std::string reads_file_a = leftover_args[0];
    const std::string reads_file_b = leftover_args[1];

    /**********************************************************************/

    /***********************************
     * STEP 1: READ IN THE DATA
     */
    std::vector<SimpleGenomicRegion> bin_boundaries;
    std::vector<double> read_bins_a;
    std::vector<double> read_bins_b;
    std::vector<double> scales;
    std::vector<size_t> reset_points;
    LoadReadsByRegion(VERBOSE, chroms_file, reads_file_a, reads_file_b,
                      deads_file, bin_size_step, bin_boundaries, read_bins_a,
                      read_bins_b, scales, reset_points, FRAGMENT_LEN,
                      BAM_FORMAT, REMOVE_JACKPOT);

    if (VERBOSE)
      std::cout << "[SELECTING BIN SIZE] ";
    if (bin_size == 0) {
      if (hideaki)
        bin_size =
          select_bin_size_hideaki(read_bins_a, scales, bin_size_step, smooth);
      else if (waterman)
        bin_size =
          select_bin_size_waterman(read_bins_a, scales, bin_size_step, smooth);
      else
        bin_size =
          select_bin_size_hideaki_emp(read_bins_a, scales, reset_points,
                                      bin_size_step, max_dead_proportion);
    }
    if (VERBOSE)
      std::cout << "bin size =  " << bin_size << '\n';

    /***********************************
     * STEP 2: BIN THE READS
     */
    AdjustBinSize(bin_boundaries, read_bins_a, read_bins_b, scales,
                  reset_points, bin_size_step, bin_size);
    RemoveDeserts(bin_boundaries, read_bins_a, read_bins_b, scales,
                  reset_points, bin_size, desert_size, max_dead_proportion);

    std::vector<double> read_bins(read_bins_a.size());
    const double max_count = bin_size;
    for (size_t i = 0; i < read_bins.size(); ++i) {
      read_bins_a[i] = std::min(read_bins_a[i], max_count);
      read_bins_b[i] = std::min(read_bins_b[i], max_count);
      read_bins[i] = read_bins_a[i] - read_bins_b[i];
    }

    std::vector<std::vector<SimpleGenomicRegion>> bin_boundaries_folded;
    expand_bins(bin_boundaries, reset_points, bin_boundaries_folded);

    // mode specific code
    if (mode == TEST_CONTROL_MODE) {
      /***********************************
       * STEP 3: ESTIMATE EMISSION PARAMS
       */

      const TwoStateScaleSplitHMM hmm(min_prob, tolerance, max_iterations,
                                      VERBOSE);
      size_t state_num = 2;
      std::vector<SplitDistro> distros;
      std::vector<std::vector<double>> trans;
      std::vector<double> start_trans, end_trans;

      if (!in_param_file.empty())
        read_param_file(in_param_file, state_num, start_trans, trans, end_trans,
                        distros);
      else {
        if (VERBOSE)
          std::cout << "[ESTIMATING PARAMETERS]" << '\n';

        fg_size = (fg_size > 0) ? fg_size : 20000;

        training_size = (training_size == 0) ? read_bins.size() : training_size;

        std::vector<double> read_bins_sample, read_bins_a_sample,
          read_bins_b_sample, scales_sample;
        std::vector<size_t> reset_points_sample;
        pick_training_sample(read_bins, read_bins_a, read_bins_b, scales,
                             reset_points, training_size, read_bins_sample,
                             read_bins_a_sample, read_bins_b_sample,
                             scales_sample, reset_points_sample);

        distros.push_back(SplitDistro(fg_name));
        distros.push_back(SplitDistro(bg_name));

        double mixing = 0;
        TwoStateSplitResolveMixture(read_bins_sample, read_bins_a_sample,
                                    read_bins_b_sample, scales_sample,
                                    MAX_INITIALIZATION_ITR, tolerance, VERBOSE,
                                    distros.front(), distros.back(), mixing);

        // void TwoStateSplitResolveMixture(
        //   const std::vector<double> &values, const std::vector<double>
        //   &vals_a, const std::vector<double> &vals_b, const
        //   std::vector<double> &scales, const std::size_t max_iterations,
        //   const double tolerance, const bool verbose, SplitDistro &fg_distro,
        //   SplitDistro &bg_distro, double &mixing) {

        set_transitions(bin_size, fg_size, mixing, start_trans, trans,
                        end_trans);

        hmm.BaumWelchTraining(read_bins_sample, read_bins_a_sample,
                              read_bins_b_sample, scales_sample,
                              reset_points_sample, start_trans, trans,
                              end_trans, distros.front(), distros.back());

        clear_training_sample(read_bins_sample, read_bins_a_sample,
                              read_bins_b_sample, scales_sample,
                              reset_points_sample);
      }

      if (!out_param_file.empty())
        write_param_file(out_param_file, state_num, start_trans, trans,
                         end_trans, distros);

      if (VERBOSE)
        report_final_values(distros, trans);

      /***********************************
       * STEP 5: DECODE THE DOMAINS
       */

      std::vector<bool> classes;
      std::vector<double> scores;
      if (USE_POSTERIOR)
        hmm.PosteriorDecoding(read_bins, scales, reset_points, start_trans,
                              trans, end_trans, distros.front(), distros.back(),
                              classes, scores);
      else
        hmm.ViterbiDecoding(read_bins, scales, reset_points, start_trans, trans,
                            end_trans, distros.front(), distros.back(),
                            classes);

      /***********************************
       * STEP 6: WRITE THE RESULTS
       */

      // make sure the output dir is valid
      output_domains(read_bins, scales, reset_points, bin_boundaries_folded,
                     hmm, distros, start_trans, trans, end_trans, classes,
                     posterior_cutoff, undef_region_cutoff, cdf_cutoff,
                     domain_file, posterior_score_file, VERBOSE);
      if (!boundary_file.empty() && boundary_file != "None")
        output_boundaries(read_bins, scales, reset_points,
                          bin_boundaries_folded, hmm, distros, start_trans,
                          trans, end_trans, classes, boundary_file,
                          boundary_score_file, VERBOSE);

      if (!read_counts_file.empty() && read_counts_file != "None") {
        write_read_counts_by_bin(bin_boundaries_folded, read_bins_a,
                                 read_bins_b, scales, classes,
                                 read_counts_file);
      }
    }
    else if (mode == TEST_TEST_MODE) {
      /***********************************
       * STEP 3: ESTIMATE EMISSION PARAMS
       */
      const ThreeStateScaleSplitHMM hmm(min_prob, tolerance, max_iterations,
                                        VERBOSE);
      size_t state_num = 3;
      std::vector<SplitDistro> distros;
      std::vector<std::vector<double>> trans;
      std::vector<double> start_trans, end_trans;

      if (!in_param_file.empty())
        read_param_file(in_param_file, state_num, start_trans, trans, end_trans,
                        distros);
      else {
        if (VERBOSE)
          std::cout << "[ESTIMATING PARAMETERS]" << '\n';

        fg_size = (fg_size > 0) ? fg_size : 6000;

        // All are using the fg name now
        training_size = (training_size == 0) ? read_bins.size() : training_size;

        std::vector<double> read_bins_sample, read_bins_a_sample,
          read_bins_b_sample, scales_sample;
        std::vector<size_t> reset_points_sample;
        pick_training_sample(read_bins, read_bins_a, read_bins_b, scales,
                             reset_points, training_size, read_bins_sample,
                             read_bins_a_sample, read_bins_b_sample,
                             scales_sample, reset_points_sample);

        distros.push_back(SplitDistro(fg_name));
        distros.push_back(SplitDistro(bg_name));
        distros.push_back(SplitDistro(fg_name));

        std::vector<double> mixing;
        ThreeStateScaleSplitResolveMixture(
          read_bins_sample, read_bins_a_sample, read_bins_b_sample,
          scales_sample, MAX_INITIALIZATION_ITR, tolerance, VERBOSE,
          distros.front(), distros[1], distros.back(), mixing);

        /***********************************
         * STEP 4: TRAIN THE HMM
         */
        set_transitions(bin_size, fg_size, mixing, start_trans, trans,
                        end_trans);

        hmm.BaumWelchTraining(
          read_bins_sample, read_bins_a_sample, read_bins_b_sample,
          scales_sample, reset_points_sample, start_trans, trans, end_trans,
          distros.front(), distros[1], distros.back());

        clear_training_sample(read_bins_sample, read_bins_a_sample,
                              read_bins_b_sample, scales_sample,
                              reset_points_sample);
      }

      if (!out_param_file.empty())
        write_param_file(out_param_file, state_num, start_trans, trans,
                         end_trans, distros);

      if (VERBOSE)
        report_final_values(distros, trans);

      /***********************************
       * STEP 5: DECODE THE DOMAINS
       */

      std::vector<size_t> classes;
      std::vector<double> scores;
      if (USE_POSTERIOR)
        hmm.PosteriorDecoding(read_bins, scales, reset_points, start_trans,
                              trans, end_trans, distros.front(), distros[1],
                              distros.back(), classes, scores);
      else
        hmm.ViterbiDecoding(read_bins, scales, reset_points, start_trans, trans,
                            end_trans, distros.front(), distros[1],
                            distros.back(), classes);

      /***********************************
       * STEP 6: WRITE THE RESULTS
       */

      // make sure the output dir is valid
      output_domains(read_bins, scales, reset_points, bin_boundaries_folded,
                     hmm, distros, start_trans, trans, end_trans, classes,
                     posterior_cutoff, undef_region_cutoff, cdf_cutoff,
                     domain_file, posterior_score_file, VERBOSE);
      if (!boundary_file.empty() && boundary_file != "None")
        output_boundaries(read_bins, scales, reset_points,
                          bin_boundaries_folded, hmm, distros, start_trans,
                          trans, end_trans, classes, boundary_file,
                          boundary_score_file, VERBOSE);

      if (!read_counts_file.empty() && read_counts_file != "None") {
        write_read_counts_by_bin(bin_boundaries_folded, read_bins_a,
                                 read_bins_b, scales, classes,
                                 read_counts_file);
      }
    }
    else {
      const auto message = R"(

Please specifc the following value for mode: if you want to compare a test
sample and a control sample and think there are two stats, use mode 2; if you
want to compare a test sample and another test sample and think thre are three
states, use mode 3

)";
      std::cerr << message;
    }
  }
  catch (std::runtime_error &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
