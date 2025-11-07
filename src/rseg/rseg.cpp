/* Copyright (C) 2011 University of Southern California
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

static constexpr auto about = R"(rseg

Segment the genome according to mapped read density.
)";

#include "Distro.hpp"
#include "EvaluateBoundaries.hpp"
#include "LoadReadsByRegion.hpp"
#include "ModelParams.hpp"
#include "ReadCounts.hpp"
#include "SelectBinSize.hpp"
#include "TwoStateScaleHMM.hpp"
#include "TwoStateScaleResolveMixture.hpp"
#include "rseg_utils.hpp"

#include "CLI11.hpp"

#include "GenomicRegion.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <print>
#include <random>

// Determines how many iterations are used during the initialization
// phase to find good starting values for the HMM
const std::size_t MAX_INITIALIZATION_ITR = 3;

static void
output_boundaries(
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const std::vector<double> &tmp_read_bins, const std::vector<double> &scales,
  const std::vector<bool> &tmp_classes,
  const std::vector<std::size_t> &reset_points, const TwoStateScaleHMM &hmm,
  const std::vector<Distro> &distros, const std::vector<double> &start_trans,
  const std::vector<std::vector<double>> &trans,
  const std::vector<double> &end_trans, const std::string &boundary_file,
  const std::string &boundary_score_file, const bool VERBOSE,
  const bool both_domain_ends) {
  static constexpr auto FDR = 0.05;

  const Distro &fg_distro = distros.front();
  const Distro &bg_distro = distros.back();

  const auto f_to_f_scores = hmm.TransitionPosteriors(
    tmp_read_bins, scales, reset_points, start_trans, trans, end_trans,
    fg_distro, bg_distro, TwoStateScaleHMM::FG_TO_FG_TRANSITION);
  const auto f_to_b_scores = hmm.TransitionPosteriors(
    tmp_read_bins, scales, reset_points, start_trans, trans, end_trans,
    fg_distro, bg_distro, TwoStateScaleHMM::FG_TO_BG_TRANSITION);
  const auto b_to_f_scores = hmm.TransitionPosteriors(
    tmp_read_bins, scales, reset_points, start_trans, trans, end_trans,
    fg_distro, bg_distro, TwoStateScaleHMM::BG_TO_FG_TRANSITION);
  const auto b_to_b_scores = hmm.TransitionPosteriors(
    tmp_read_bins, scales, reset_points, start_trans, trans, end_trans,
    fg_distro, bg_distro, TwoStateScaleHMM::BG_TO_BG_TRANSITION);

  std::vector<double> tmp_boundary_scores(std::size(f_to_b_scores));
  std::vector<std::uint32_t> transitions(1, 0);
  for (auto i = 1u; i < std::size(tmp_classes); ++i)
    if (tmp_classes[i] != tmp_classes[i - 1])
      transitions.push_back(i);
  transitions.push_back(std::size(tmp_classes));

  // ADS: check that number of transitions is never more than number of bins
  auto j = 0u;
  for (auto i = 0u; i < std::size(f_to_b_scores); ++i) {
    if (i > (transitions[j] + transitions[j + 1]) / 2u)
      ++j;
    tmp_boundary_scores[i] =
      tmp_classes[transitions[j]] ? b_to_f_scores[i] : f_to_b_scores[i];
  }

  std::vector<std::vector<double>> boundary_scores;
  expand_bins(tmp_boundary_scores, reset_points, boundary_scores);

  // generate control sample
  std::random_device rd{};
  std::mt19937 g{rd()};

  std::vector<double> read_counts_control(tmp_read_bins);
  for (std::size_t i = 0; i < std::size(reset_points) - 1; ++i)
    std::shuffle(std::begin(read_counts_control) + reset_points[i],
                 std::begin(read_counts_control) + reset_points[i + 1], g);

  const auto b_to_f_scores_control = hmm.TransitionPosteriors(
    read_counts_control, scales, reset_points, start_trans, trans, end_trans,
    fg_distro, bg_distro, TwoStateScaleHMM::BG_TO_FG_TRANSITION);

  const auto f_to_b_scores_control = hmm.TransitionPosteriors(
    read_counts_control, scales, reset_points, start_trans, trans, end_trans,
    fg_distro, bg_distro, TwoStateScaleHMM::FG_TO_BG_TRANSITION);

  std::vector<double> boundary_scores_control(std::size(f_to_b_scores_control),
                                              0);

  for (std::size_t i = 0; i < std::size(boundary_scores_control); ++i)
    boundary_scores_control[i] =
      std::max(f_to_b_scores_control[i], b_to_f_scores_control[i]);

  std::sort(std::begin(boundary_scores_control),
            std::end(boundary_scores_control));
  const std::size_t bc_idx = std::size(boundary_scores_control) * (1.0 - FDR);
  const double cutoff = boundary_scores_control[bc_idx];

  // read_counts_control.clear();
  // read_counts_control.shrink_to_fit();
  // b_to_f_scores_control.clear();
  // b_to_f_scores_control.shrink_to_fit();
  // f_to_b_scores_control.clear();
  // f_to_b_scores_control.shrink_to_fit();
  // boundary_scores_control.clear();
  // boundary_scores_control.shrink_to_fit();

  // finished generating control sample

  std::vector<GenomicRegion> boundaries;
  BoundEval be(1, 1);
  if (both_domain_ends)
    be.evaluate(bin_bounds, reset_points, tmp_classes, tmp_boundary_scores,
                f_to_f_scores, f_to_b_scores, b_to_f_scores, b_to_b_scores,
                cutoff, both_domain_ends, boundaries);
  else
    be.evaluate(bin_bounds, reset_points, /*tmp_classes,*/ tmp_boundary_scores,
                f_to_f_scores, f_to_b_scores, b_to_f_scores, b_to_b_scores,
                cutoff, boundaries);

  // write result files
  if (VERBOSE)
    std::cout << "Boundary file: " + boundary_file << '\n';
  WriteBEDFile(boundary_file, boundaries);

  if (!boundary_score_file.empty()) {
    if (VERBOSE)
      std::cout << "Boundary score file: " + boundary_score_file << '\n';
    write_wigfile(boundary_scores, bin_bounds, boundary_score_file);
  }
}

static void
output_domains(const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
               const std::vector<double> &tmp_read_bins,
               const std::vector<double> &tmp_scales,
               const std::vector<bool> &tmp_classes,
               std::vector<double> &tmp_scores,
               const std::vector<std::size_t> &reset_points,
               const TwoStateScaleHMM &hmm, const std::vector<Distro> &distros,
               const std::vector<double> &start_trans,
               const std::vector<std::vector<double>> &trans,
               const std::vector<double> &end_trans,
               const double posterior_cutoff,
               const std::size_t undef_region_cutoff, const double cdf_cutoff,
               const std::string &output_file,
               const std::string &posterior_score_file, const bool VERBOSE) {

  // Obtain the scores for the current domain class
  if (tmp_scores.size() == 0)
    tmp_scores = hmm.PosteriorScores(
      tmp_read_bins, tmp_scales, reset_points, start_trans, trans, end_trans,
      distros.front(), distros.back(), tmp_classes);

  std::vector<std::vector<double>> read_bins, scores, scales;
  std::vector<std::vector<bool>> classes;
  expand_bins(tmp_read_bins, reset_points, read_bins);
  expand_bins(tmp_scores, reset_points, scores);
  expand_bins(tmp_scales, reset_points, scales);
  expand_bins(tmp_classes, reset_points, classes);

  std::vector<std::vector<GenomicRegion>> domains;
  build_domains(bin_bounds, classes, scores, posterior_cutoff, domains,
                undef_region_cutoff);

  // output domains
  pick_domains(bin_bounds, read_bins, scales, distros, domains, cdf_cutoff);

  write_bed_file(domains, output_file);
  if (VERBOSE)
    std::cout << "Domains file: " + output_file << '\n';

  if (!posterior_score_file.empty()) {
    std::size_t k = 0;
    for (std::size_t i = 0; i < scores.size(); ++i)
      for (std::size_t j = 0; j < scores[i].size(); ++j) {
        scores[i][j] == tmp_classes[k] ? scores[i][j] : 1 - scores[i][j];
        ++k;
      }
    write_wigfile(scores, bin_bounds, posterior_score_file);
    if (VERBOSE)
      std::cout << "Bin score file: " + posterior_score_file << '\n';
  }
}

int
main(int argc, char *argv[]) {
  static constexpr auto usage = "Usage: rseg [options]";
  // names of emission distributions to use
  static constexpr auto fg_name = "nbd";
  static constexpr auto bg_name = "nbd";
  static constexpr auto both_domain_ends{true};

  try {
    // file names
    std::string reads_file;
    std::string deads_file;
    std::string chroms_file;
    std::string in_param_file;
    std::string out_param_file;

    // expected size of a domain
    double fg_size = 20000;

    // flags
    bool USE_POSTERIOR = false;
    bool REMOVE_JACKPOT{false};
    bool VERBOSE = false;
    bool BAM_FORMAT = false;

    std::string output_file;
    std::string posterior_score_file;
    std::string boundary_file;
    std::string boundary_score_file;
    std::string read_counts_file;

    std::size_t desert_size = 20000;
    std::size_t bin_size_step = 50;
    std::size_t bin_size = 0;
    std::size_t FRAGMENT_LEN = 0;
    bool waterman = false;
    bool hideaki = false;
    bool hideaki_emp = false;
    bool smooth = true;
    std::size_t max_iterations = 20;
    double tolerance = 1e-20;
    double min_prob = 1e-20;

    // posterior theshhold above which a bin is considerd belonging to a state
    double posterior_cutoff = 0.95;

    // if an undefined region larger then this value, leave it as is
    std::size_t undef_region_cutoff = 3000;

    double cdf_cutoff = 0.1;

    double max_dead_proportion = 0.5;

    CLI::App app{about};
    argv = app.ensure_utf8(argv);
    app.usage(usage);

    // clang-format off
    app.set_help_flag("-h,--help", "Print a detailed help message and exit");
    app.add_option("-o,--out", output_file, "output file")
      ->required();
    app.add_option("--score", posterior_score_file, "Posterior scores file");
    app.add_option("--readcount", read_counts_file, "readcounts file");
    app.add_option("--boundary", boundary_file, "domain boundary file");
    app.add_option("--boundary-score", boundary_score_file, "boundary transition scores file");
    app.add_option("-r,--reads", reads_file, "mapped reads file (BED or BAM format)")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-c,--chrom", chroms_file, "file with chromosome sizes (BED format)")
      ->required()
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-d,--deadzones", deads_file, "file of deadzones (BED format)")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_flag("--bam", BAM_FORMAT, "Input reads file is BAM format");
    app.add_option("--param-in", in_param_file, "Input parameters file")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("--param-out", out_param_file, "Output parameters file");
    app.add_option("-i,--maxitr", max_iterations, "maximum iterations for training");
    app.add_option("-b,--bin-size", bin_size, "bin size (default: based on data)");
    app.add_option("--bin-step", bin_size_step, "minimum bin size");
    app.add_flag("--no-jackpots", REMOVE_JACKPOT, "remove jackpot reads");
    app.add_option("--fragment_length", FRAGMENT_LEN,
                   "Extend reads to fragment length (default not to extend)");
    app.add_flag("--Waterman", waterman, "use Waterman's method for bin size");
    app.add_flag("--Hideaki", hideaki, "use Hideaki's method for bin size");
    app.add_flag("--Hideaki-emp", hideaki_emp, "use Hideaki's empirical method (default)");
    app.add_flag("--smooth", smooth, "Indicate whether the rate curve is assumed smooth");
    app.add_option("--max-dead", max_dead_proportion, "max deadzone proportion for retained bins");
    app.add_option("-s,--domain-size", fg_size, "expected domain size");
    app.add_option("-S,--desert", desert_size, "desert size");
    // app.add_option("-F,--fg", fg_name, "foreground emission distribution");
    // app.add_option("-B,--bg", bg_name, "background emission distribution");
    app.add_flag("-P,--posterior", USE_POSTERIOR, "use posterior decoding (default: Viterbi)");
    app.add_option("--posterior-cutoff", posterior_cutoff, "posterior cutoff significance");
    app.add_option("--undefined", undef_region_cutoff, "min size of unmappable region");
    app.add_option("--cutoff", cdf_cutoff, "cutoff in cdf for identified domains");
    app.add_flag("-v,--verbose", VERBOSE, "print more info");
    // clang-format on

    if (argc == 1) {
      std::println("{}", app.help());
      return EXIT_FAILURE;
    }

    CLI11_PARSE(app, argc, argv);

    if (VERBOSE)
      std::cout << "[PROCESSING] " << strip_path(reads_file) << '\n';

    // read in the data
    std::vector<SimpleGenomicRegion> bin_boundaries;
    std::vector<double> read_bins;
    std::vector<double> scales;
    std::vector<std::size_t> reset_points;
    LoadReadsByRegion(VERBOSE, chroms_file, reads_file, deads_file,
                      bin_size_step, bin_boundaries, read_bins, scales,
                      reset_points, FRAGMENT_LEN, BAM_FORMAT, REMOVE_JACKPOT);

    if (VERBOSE)
      std::cout << "[SELECTING BIN SIZE] ";
    bin_size = [&]() {
      if (bin_size != 0)
        return bin_size;
      if (hideaki)
        return select_bin_size_hideaki(read_bins, scales, bin_size_step,
                                       smooth);
      if (waterman)
        return select_bin_size_waterman(read_bins, scales, bin_size_step,
                                        smooth);
      return select_bin_size_hideaki_emp(read_bins, scales, reset_points,
                                         bin_size_step, max_dead_proportion);
    }();
    if (VERBOSE)
      std::cout << "bin size =  " << bin_size << '\n';

    if (!std::all_of(std::cbegin(scales), std::cend(scales),
                     [](const auto x) { return x >= 0.0 && x <= 1.0; }))
      throw std::runtime_error("not all scales correct");

    /// make bins of reads
    AdjustBinSize(bin_boundaries, read_bins, scales, reset_points,
                  bin_size_step, bin_size);

    if (!std::all_of(std::cbegin(scales), std::cend(scales),
                     [](const auto x) { return x >= 0.0 && x <= 1.0; }))
      throw std::runtime_error("not all scales correct after adjustment");

    RemoveDeserts(bin_boundaries, read_bins, scales, reset_points, bin_size,
                  desert_size, max_dead_proportion);

    const double max_count = bin_size;
    for (std::size_t i = 0; i < read_bins.size(); ++i)
      read_bins[i] = std::min(read_bins[i], max_count);

    if (!std::all_of(std::cbegin(read_bins), std::cend(read_bins),
                     [](const auto x) { return x >= 0.0; }))
      throw std::runtime_error("non-negative bins found");

    std::vector<std::vector<SimpleGenomicRegion>> bin_boundaries_folded;
    expand_bins(bin_boundaries, reset_points, bin_boundaries_folded);

    // estimate emission params
    const TwoStateScaleHMM hmm(min_prob, tolerance, max_iterations, VERBOSE);
    std::size_t state_num = 2;
    std::vector<Distro> distros;
    std::vector<std::vector<double>> trans;
    std::vector<double> start_trans, end_trans;

    if (!in_param_file.empty())
      read_param_file(in_param_file, state_num, start_trans, trans, end_trans,
                      distros);
    else {
      if (VERBOSE)
        std::cout << "[ESTIMATING PARAMETERS]" << '\n';

      distros.push_back(Distro(fg_name));
      distros.push_back(Distro(bg_name));

      assert(read_bins.size() == scales.size());

      double mixing = 0;
      TwoStateResolveMixture(read_bins, scales, MAX_INITIALIZATION_ITR,
                             tolerance, VERBOSE, distros.front(),
                             distros.back(), mixing);

      set_transitions(bin_size, fg_size, mixing, start_trans, trans, end_trans);

      hmm.BaumWelchTraining(read_bins, scales, reset_points, start_trans, trans,
                            end_trans, distros.front(), distros.back());
    }

    if (!out_param_file.empty())
      write_param_file(out_param_file, state_num, start_trans, trans, end_trans,
                       distros);

    if (VERBOSE)
      report_final_values(distros, trans);

    // decode the domains
    std::vector<bool> classes;
    std::vector<double> scores;
    if (USE_POSTERIOR)
      hmm.PosteriorDecoding(read_bins, scales, reset_points, start_trans, trans,
                            end_trans, distros.front(), distros.back(), classes,
                            scores);
    else
      hmm.ViterbiDecoding(read_bins, scales, reset_points, start_trans, trans,
                          end_trans, distros.front(), distros.back(), classes);

    // write the results
    // make sure the output dir is valid
    output_domains(bin_boundaries_folded, read_bins, scales, classes, scores,
                   reset_points, hmm, distros, start_trans, trans, end_trans,
                   posterior_cutoff, undef_region_cutoff, cdf_cutoff,
                   output_file, posterior_score_file, VERBOSE);
    if (!boundary_file.empty())
      output_boundaries(bin_boundaries_folded, read_bins, scales, classes,
                        reset_points, hmm, distros, start_trans, trans,
                        end_trans, boundary_file, boundary_score_file, VERBOSE,
                        both_domain_ends);

    if (!read_counts_file.empty())
      write_read_counts_by_bin(bin_boundaries_folded, read_bins, scales,
                               classes, read_counts_file);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
