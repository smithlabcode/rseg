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

#include <numeric>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <iostream>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "Distro.hpp"
#include "ReadCounts.hpp"
#include "TwoStateResolveMixture.hpp"
#include "TwoStateScaleHMM.hpp"
#include "EvaluateBoundaries.hpp"
#include "LoadReadsByRegion.hpp"
#include "SelectBinSize.hpp"
#include "OptionParser.hpp"
#include "rseg_utils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;


// Determines how many iterations are used during the initialization
// phase to find good starting values for the HMM
const size_t MAX_INITIALIZATION_ITR = 3;

static void
output_boundaries(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
                  const vector<double> &tmp_read_bins,
                  const vector<double> &scales,
                  const vector<bool> &tmp_classes,
                  const vector<size_t> &reset_points,
                  const TwoStateScaleHMM &hmm,
                  const vector<Distro> &distros,
                  const vector<double> &start_trans,
                  const vector<vector<double> > &trans,
                  const vector<double> &end_trans,
                  const string &dataset_name, const string &outdir, 
                  const bool VERBOSE, const bool WRITE_TRACKS,
                  const bool Both_Domain_Ends = true) {

  static const double FDR = 0.05;
  
  static const string BED_SUFF = string(".bed");
  static const string WIG_SUFF = string(".wig");
  static const string DOMAINS_TAG = string("-domains");
  static const string BOUNDARY_TAG = string("-boundaries");
  static const string SCORES_TAG = string("-scores");

  const Distro &fg_distro = distros.front();
  const Distro &bg_distro = distros.back();
  vector<double> f_to_f_scores, f_to_b_scores, b_to_f_scores, b_to_b_scores;
  hmm.TransitionPosteriors(tmp_read_bins, scales, reset_points,
			   start_trans, trans, end_trans, fg_distro, bg_distro,
			   TwoStateScaleHMM::FG_TO_FG_TRANSITION, f_to_f_scores);
  hmm.TransitionPosteriors(tmp_read_bins, scales, reset_points,
			   start_trans, trans, end_trans, fg_distro, bg_distro,
			   TwoStateScaleHMM::FG_TO_BG_TRANSITION, f_to_b_scores);
  hmm.TransitionPosteriors(tmp_read_bins, scales, reset_points,
			   start_trans, trans, end_trans, fg_distro, bg_distro,
			   TwoStateScaleHMM::BG_TO_FG_TRANSITION, b_to_f_scores);
  hmm.TransitionPosteriors(tmp_read_bins, scales, reset_points,
			   start_trans, trans, end_trans, fg_distro, bg_distro,
			   TwoStateScaleHMM::BG_TO_BG_TRANSITION, b_to_b_scores);

    
  vector<double> tmp_boundary_scores(f_to_b_scores.size());
  vector<int> transitions;
  for (size_t i = 0; i < tmp_classes.size(); ++i)
    if (i == 0 || tmp_classes[i] != tmp_classes[i - 1])
      transitions.push_back(i);
    
  transitions.push_back(tmp_classes.size());
  size_t j = 0;
  for (int i = 0; static_cast<size_t>(i) < f_to_b_scores.size(); ++i) {
    if (abs(i - transitions[j]) > abs(i - transitions[j + 1]))
      ++j;
    tmp_boundary_scores[i] = (tmp_classes[transitions[j]]) ?
      b_to_f_scores[i] : f_to_b_scores[i];
  }
  
  vector<vector<double> > boundary_scores;
  expand_bins(tmp_boundary_scores, reset_points, boundary_scores);
    

  //// generate control sample
  vector<double> read_counts_control(tmp_read_bins);
  for (size_t i = 0; i < reset_points.size() - 1; ++i)
    std::random_shuffle(read_counts_control.begin() + reset_points[i],
			read_counts_control.begin() + reset_points[i + 1]);

  vector<double> b_to_f_scores_control, f_to_b_scores_control;
  hmm.TransitionPosteriors(read_counts_control, scales, reset_points,
			   start_trans, trans, end_trans, fg_distro, bg_distro,
			   TwoStateScaleHMM::BG_TO_FG_TRANSITION, b_to_f_scores_control);
  hmm.TransitionPosteriors(read_counts_control,  scales, reset_points,
			   start_trans, trans, end_trans, fg_distro, bg_distro, 
			   TwoStateScaleHMM::FG_TO_BG_TRANSITION, f_to_b_scores_control);

  vector<double> boundary_scores_control(f_to_b_scores_control.size(), 0);

  for (size_t i = 0; i < boundary_scores_control.size(); ++i)
    boundary_scores_control[i] = std::max(f_to_b_scores_control[i],
					  b_to_f_scores_control[i]);

  std::sort(boundary_scores_control.begin(), boundary_scores_control.end());
  const size_t bc_idx = static_cast<size_t>(boundary_scores_control.size()*(1 - FDR));
  const double cutoff = boundary_scores_control[bc_idx];
  
  read_counts_control.clear();
  b_to_f_scores_control.clear();
  f_to_b_scores_control.clear();
  boundary_scores_control.clear();
  //// finish generating control sample

  vector<GenomicRegion> boundaries;
  BoundEval be(1, 1);
  if (Both_Domain_Ends)
    be.evaluate(bin_bounds, reset_points, tmp_classes, tmp_boundary_scores,
		f_to_f_scores, f_to_b_scores, b_to_f_scores, b_to_b_scores,
		cutoff, Both_Domain_Ends, boundaries);
  else be.evaluate(bin_bounds, reset_points, tmp_classes, tmp_boundary_scores,
		   f_to_f_scores, f_to_b_scores, b_to_f_scores, b_to_b_scores,
		   cutoff, boundaries);
  
  // write result files
  const string boundary_filename(path_join(outdir, dataset_name + 
					   BOUNDARY_TAG + BED_SUFF));
  if (VERBOSE)
    cerr << "Boundary file: " + boundary_filename << endl;
  WriteBEDFile(boundary_filename, boundaries);
  
  if (WRITE_TRACKS) {
    const string bound_scores_filename(path_join(outdir, dataset_name + 
						 "-boundary-scores" + WIG_SUFF));
    if (VERBOSE)
      cerr << "Boundary score file: " + bound_scores_filename << endl;
    write_wigfile(boundary_scores, bin_bounds, bound_scores_filename);
  }
}


static void
output_domains(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
	       const vector<double> &tmp_read_bins,
               const vector<double> &tmp_scales,
	       const vector<bool> &tmp_classes,
               vector<double> &tmp_scores,
	       const vector<size_t> &reset_points,
	       const TwoStateScaleHMM &hmm,
	       const vector<Distro> &distros,
	       const vector<double> &start_trans,
	       const vector<vector<double> > &trans,
               const vector<double> &end_trans,
               const double posterior_cutoff,
               const size_t undef_region_cutoff,
               const double cdf_cutoff, 
	       const string dataset_name, const string outdir, 
	       const bool VERBOSE, const bool WRITE_TRACKS) {
  
  static const string BED_SUFF = string(".bed");
  static const string WIG_SUFF = string(".wig");
  static const string DOMAINS_TAG = string("-domains");
  static const string BOUNDARY_TAG = string("-boundaries");
  static const string SCORES_TAG = string("-scores");
  
  // Obtain the scores for the current domain class
  if (tmp_scores.size() == 0)
    hmm.PosteriorScores(tmp_read_bins, tmp_scales, reset_points,
			start_trans, trans, end_trans, 
			distros.front(), distros.back(), tmp_classes, tmp_scores);
  
  vector<vector<double> > read_bins, scores, scales;
  vector<vector<bool> > classes;
  expand_bins(tmp_read_bins, reset_points, read_bins);
  expand_bins(tmp_scores, reset_points, scores);
  expand_bins(tmp_scales, reset_points, scales);
  expand_bins(tmp_classes, reset_points, classes);
  
  vector<vector<GenomicRegion> > domains;
  build_domains(bin_bounds, classes, scores, posterior_cutoff,
		domains, undef_region_cutoff);
    
  // output domains
  pick_domains(bin_bounds, read_bins, scales, distros,
	       domains, cdf_cutoff);
    
  const string domain_file_name =
    path_join(outdir ,dataset_name + DOMAINS_TAG + BED_SUFF);
  write_bed_file(domains, domain_file_name);
  if (VERBOSE)
    cerr << "Domains file: " + domain_file_name << endl;
    
  if (WRITE_TRACKS)
    {
      const string scores_file_name(path_join(outdir, dataset_name + 
					      SCORES_TAG + WIG_SUFF));
      write_wigfile(scores, bin_bounds, scores_file_name);
      if (VERBOSE)
	cerr << "Bin score file: " + scores_file_name << endl;
    }
}


int
main(int argc, const char **argv) {
  try {
    // file names
    string reads_file;
    string deads_file;
    string chroms_file;

    string outdir = ".";
    string tmp_dataset_name;
  
    // expected size of a domain
    double fg_size = 20000;
        
    // flags
    bool use_posterior = true;
    bool REMOVE_JACKPOT = true;
    bool VERBOSE = false;
    bool WRITE_BOUNDARY = false;
    bool WRITE_TRACKS = false;
    bool PRINT_READCOUNTS = false;
    bool Both_Domain_Ends = true;
    
    // names of statistical distributions to use
    string fg_name = "nbd";
    string bg_name = "nbd";
  
    size_t desert_size = 20000;
    size_t bin_size_step = 100;
    size_t bin_size = 0;
    bool waterman = false;
    bool hideaki = false;
    bool hideaki_emp = false;
    bool smooth = true;
    size_t max_iterations = 20;
    double tolerance = 1e-20;
    double min_prob = 1e-20;
    
    // the posterior theshhold above which a bin is considerd belonging to a state
    double posterior_cutoff = 0.95; 
    
    // if an undefined region larger then this value, leave it as is
    size_t undef_region_cutoff = 3000; 
        
    double cdf_cutoff =  0.1;
    
    double max_dead_proportion = 0.5;
    
    ////////////////////// PARSING COMMAND LINE OPTIONS /////////////////////////
    OptionParser opt_parse(basename(argv[0]),
                           "This program segments genome according "
                           "to mapped read density", "BED_file");
    opt_parse.add_opt("output-dir", 'o', "Output directory name (default CWD)", 
		      false, outdir);
    opt_parse.add_opt("boundary", '\0', "Write boundary file", 
		      false, WRITE_BOUNDARY);
    opt_parse.add_opt("tracks", '\0', "Whether write additional browser tracks", 
		      false, WRITE_TRACKS);
    opt_parse.add_opt("read-counts", '\0', "Write reads counts "
		      "file in each bin", false, PRINT_READCOUNTS);
    opt_parse.add_opt("name", '\0', "Name of dataset (default: filename)", 
		      false, tmp_dataset_name);
    opt_parse.add_opt("chrom", 'c', "Name of the file with sizes of chromosomes", 
		      true, chroms_file);
    opt_parse.add_opt("deadzone-file", 'd', "Filename of deadzones", false, deads_file);
    opt_parse.add_opt("domain-size", 's', "Expected size of domain (Default 20000)", 
		      false, fg_size);
    opt_parse.add_opt("bin-size", 'b', "Size of bins (default depends on # of reads)", 
		      false, bin_size);
    opt_parse.add_opt("bin-size-step", '\0',
              "Intial bin size when reading in raw reads (default 100)", 
		      false, bin_size_step);
    opt_parse.add_opt("not-remove-jackpot", '\0', "Do not remove duplicate reads", 
		      false, REMOVE_JACKPOT);
    opt_parse.add_opt("Waterman", '\0', "using Waterman's method for bin size", 
		      false, waterman);
    opt_parse.add_opt("Hideaki", '\0', "Using Hideaki's method for bin size", 
		      false, hideaki);
    opt_parse.add_opt("Hideaki-emp", '\0', "Using Hideaki's empirical method for "
		      "bin size (default)", false, hideaki_emp);
    opt_parse.add_opt("smooth", '\0', "Indicate whether the rate curve is assumed smooth", 
		      false, smooth);
    opt_parse.add_opt("max-deadzone-prop", '\0',
		      "Maximum deadzone proportion allowed for retened bins",
		      false, max_dead_proportion);
    opt_parse.add_opt("desert-size", 'S', "Desert size", false, desert_size);
    opt_parse.add_opt("iteration", 'i', "Maximum number of iterations for "
		      "HMM training", false, max_iterations);
    opt_parse.add_opt("fg", 'F', "foreground emission distribution name", false, fg_name);
    opt_parse.add_opt("bg", 'B', "background emission distribution name", false, bg_name);
    opt_parse.add_opt("posterior", 'P', "Options for posterior decoding "
		      "(default Viterbi)", false, use_posterior);
    opt_parse.add_opt("posterior-cutoff", '\0', "Posterior threshold for "
		      "signigicant bins", false, posterior_cutoff);
    opt_parse.add_opt("undef-region-cutoff", '\0', "Minimum size of "
		      "undefined region", false, undef_region_cutoff);
    opt_parse.add_opt("cdf-cutoff", '\0', "Cutoff of cumulative probability "
		      "for a true fg domain", false, cdf_cutoff); 
    // opt_parse.add_opt("tolerance", '\0', "Tolerance for convergence", 
	// 	      false, tolerance);
    // opt_parse.add_opt("min_prob", '\0', "Minimum probability value", 
	// 	      false, min_prob);
    opt_parse.add_opt("verbose", 'v', "Print more running information", 
		      false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
	
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    reads_file = leftover_args.front();
    if (reads_file.size() == 0) {
      cerr << "ERROR: input file name required" << endl;
      return EXIT_FAILURE;
    }
    if (chroms_file.size() == 0) {
      cerr << "ERROR: chromsome sizes required (BED format file)" << endl;
      return EXIT_FAILURE;
    }
    const string dataset_name = (tmp_dataset_name.size()) ? 
      tmp_dataset_name : strip_path_and_bed_suffix(reads_file.c_str());
    
    if (VERBOSE)
      cerr << "[PROCESSING] " <<  dataset_name << endl;
    
    /***********************************
     * STEP 1: READ IN THE DATA
     */

    vector<SimpleGenomicRegion> bin_boundaries;
    vector<double> read_bins;
    vector<double> scales;
    vector<size_t> reset_points;
    LoadReadsByRegion(VERBOSE, chroms_file, reads_file, deads_file, 
                      bin_size_step, bin_boundaries, read_bins,
                      scales, reset_points, REMOVE_JACKPOT);

    if (VERBOSE)
        cerr << "[SELECTING BIN SIZE] ";
    if (bin_size == 0) {
        if (hideaki)
            bin_size = select_bin_size_hideaki(
                read_bins, scales, bin_size_step, smooth);
        else if (waterman)
            bin_size = select_bin_size_waterman(
                read_bins, scales, bin_size_step, smooth);
        else 
            bin_size = select_bin_size_hideaki_emp(
                read_bins, scales, reset_points, bin_size_step,
                max_dead_proportion);
    }
    if (VERBOSE) cerr << "bin size =  " << bin_size << endl;
    
    /***********************************
     * STEP 2: BIN THE READS
     */ 
    AdjustBinSize(bin_boundaries, read_bins, scales, reset_points,
                  bin_size_step, bin_size);
    RemoveDeserts(bin_boundaries, read_bins, scales, reset_points,
                  bin_size, desert_size, max_dead_proportion);

    const double max_count = bin_size;
    for (size_t i = 0; i < read_bins.size(); ++i)
      read_bins[i] = min(read_bins[i], max_count);
    
    vector<vector<SimpleGenomicRegion> > bin_boundaries_folded;
    expand_bins(bin_boundaries, reset_points, bin_boundaries_folded);

    /***********************************
     * STEP 3: ESTIMATE EMISSION PARAMS
     */ 
    if (VERBOSE)
      cerr << "[ESTIMATIN PARAMETERS]" << endl;
    
    vector<Distro> distros;
    distros.push_back(Distro(fg_name));
    distros.push_back(Distro(bg_name));

    assert(read_bins.size() == scales.size());

    double mixing = 0;
    TwoStateResolveMixture(read_bins, scales, MAX_INITIALIZATION_ITR, 
			   tolerance, VERBOSE,
			   distros.front(), distros.back(), mixing);
    
    /***********************************
     * STEP 4: TRAIN THE HMM
     */
    
    vector<vector<double> > trans;
    vector<double> start_trans, end_trans;
    set_transitions(bin_size, fg_size, mixing, VERBOSE,
		    start_trans, trans, end_trans);
    
    const TwoStateScaleHMM hmm(min_prob, tolerance, max_iterations, VERBOSE);
    hmm.BaumWelchTraining(read_bins, scales, reset_points, start_trans, trans, 
			  end_trans, distros.front(), distros.back());
    
    if (VERBOSE)
      report_final_values(distros, start_trans, trans, end_trans);
    
    /***********************************
     * STEP 5: DECODE THE DOMAINS
     */
    vector<bool> classes;
    vector<double> scores;
    if (use_posterior)
      hmm.PosteriorDecoding(read_bins, scales, reset_points,
			    start_trans, trans, end_trans,
			    distros.front(), distros.back(), classes, scores);
    else
      hmm.ViterbiDecoding(read_bins, scales, reset_points,
			  start_trans, trans, end_trans,
			  distros.front(), distros.back(), classes);
    
    /***********************************
     * STEP 6: WRITE THE RESULTS
     */
    // make sure the output dir is valid
    chk_and_mk_dirs(outdir);
    
    output_domains(bin_boundaries_folded,
		   read_bins, scales, classes, scores, reset_points,
		   hmm, distros, start_trans, trans, end_trans, 
		   posterior_cutoff, undef_region_cutoff, cdf_cutoff,
		   dataset_name, outdir.c_str(), VERBOSE, WRITE_TRACKS);
    if (WRITE_BOUNDARY)
        output_boundaries(bin_boundaries_folded,
              read_bins, scales, classes, reset_points,
		      hmm, distros, start_trans, trans, end_trans, 
		      dataset_name, outdir,
		      VERBOSE, WRITE_TRACKS, Both_Domain_Ends);
    
    if (PRINT_READCOUNTS) {
      const string read_counts_file_name = 
	path_join(outdir, dataset_name + "-counts.bed");
      write_read_counts_by_bin(bin_boundaries_folded, read_bins, scales,
                               classes, read_counts_file_name, VERBOSE);
    }
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
