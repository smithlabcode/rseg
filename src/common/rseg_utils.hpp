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

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <utility>

#include <cmath>
#include <cassert>

#include "GenomicRegion.hpp"
#include "Distro.hpp"
#include "SplitDistro.hpp"


using std::vector;
using std::string;
using std::cin;
using std::cerr;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;

double 
get_mean(const SplitDistro &distro);

void
pick_training_sample(const vector<double> &read_bins,
                     const vector<double> &read_bins_a,
                     const vector<double> &read_bins_b,
                     const vector<double> &scales,
                     const vector<size_t> &reset_points,
                     const size_t training_sample_size,
                     vector<double> &read_bins_sample,
                     vector<double> &read_bins_a_sample,
                     vector<double> &read_bins_b_sample,
                     vector<double> &scales_sample,
                     vector<size_t> &reset_points_sample);
void
clear_training_sample(vector<double> &read_bins_sample,
		      vector<double> &read_bins_a_sample,
		      vector<double> &read_bins_b_sample,
		      vector<double> &scales_sample,
                      vector<size_t> &reset_points_sample);

double
get_max_count_cutoff(const vector<double> &values);

void
benjamini_hochberg_procedure(const vector<double> &scores,
			     const vector<double> &scores_bg,
			     vector< pair<double, size_t> > &p_values,
			     double &cutoff,
			     const double fdr);

size_t
get_bin_size(const vector<SimpleGenomicRegion> &regions,
	     const vector<vector<SimpleGenomicRegion> > &reads,
	     const vector<vector<SimpleGenomicRegion> > &deads, int VERBOSE);

void
set_transitions(const size_t bin_size, const double fg_size,
		const vector<double> &mixing,
		const bool VERBOSE,
		vector<double> &start_trans, 
		vector<vector<double> > &trans,
		vector<double> &end_trans);
void
set_transitions(const size_t bin_size, const double fg_size,
		const double mixing, const bool VERBOSE,
		vector<double> &start_trans, 
		vector<vector<double> > &trans,
		vector<double> &end_trans);

void
report_final_values(const vector<Distro> &distros,
		    const vector<double> &start_trans,
		    const vector<vector<double> > &trans,
		    const vector<double> &end_trans);
void
report_final_values(const vector<SplitDistro> &distros,
		    const vector<double> &start_trans,
		    const vector<vector<double> > &trans,
		    const vector<double> &end_trans);

void
remove_duplicate_reads(vector<SimpleGenomicRegion>  &reads);

void
remove_duplicate_reads(vector< vector<SimpleGenomicRegion> >  &reads);

void
chk_and_mk_dirs(const string & path);

void
write_read_counts_by_bin(const vector< vector<SimpleGenomicRegion> > &bin_boundaries,
                         const vector<double> &read_bins,
                         const vector<bool> &classes,
                         const string &file_name,
                         const bool VERBOSE = false);

void
write_read_counts_by_bin(const vector< vector<SimpleGenomicRegion> > &bin_boundaries,
                         const vector<double> &read_bins,
                         const vector<double> &scales,
                         const vector<bool> &classes,
                         const string &file_name,
                         const bool VERBOSE = false);
void
write_read_counts_by_bin(const vector< vector<SimpleGenomicRegion> > &bin_boundaries,
                         const vector<double> &read_bins,
                         const vector<double> &read_bins_a,
                         const vector<double> &read_bins_b,
                         const vector<bool> &classes,
                         const string &file_name,
                         const bool VERBOSE = false);
void
write_read_counts_by_bin(const vector< vector<SimpleGenomicRegion> > &bin_boundaries,
                         const vector<double> &read_bins,
                         const vector<double> &read_bins_a,
                         const vector<double> &read_bins_b,
                         const vector<size_t> &classes,
                         const string &file_name,
                         const bool VERBOSE = false);

void
elim_empty_regions(vector<SimpleGenomicRegion> &regions, 
		   vector<vector<SimpleGenomicRegion> > &bounds, 
		   vector<vector<double> > &read_bins);

void
elim_empty_regions(vector<SimpleGenomicRegion> &regions, 
		   vector<vector<SimpleGenomicRegion> > &bounds, 
		   vector<vector<double> > &read_bins,
                   vector< vector<double> > &scales);

void
elim_empty_regions(vector<SimpleGenomicRegion> &regions, 
		   vector<vector<SimpleGenomicRegion> > &bounds, 
		   vector<vector<double> > &read_bins_a, 
		   vector<vector<double> > &read_bins_b,
                   vector<vector<double> > &scales);

string
strip_path_and_bed_suffix(const string &full_path);

void
write_wigfile(const vector<vector<double> > &scores,
	      const vector<vector<SimpleGenomicRegion> > &bin_bounds,
	      const string &wigfile_name);
void
write_bed_file(const vector<vector<GenomicRegion> > &regions,
	       const string &bed_file);

template <class T> void
expand_bins(const vector<T> &tmp_bins, const vector<size_t> &resets,
	    vector<vector<T> > &bins);

void
collapse_read_bins(const vector<vector<double> > &tmp_read_bins, 
		   vector<double> &read_bins,
		   vector<size_t> &reset_points);
void
build_domains(const vector<vector<SimpleGenomicRegion> > &bins,
              const vector<vector<bool> > &classes,
              const vector<vector<double> > &scores, // posterior score of classes[i]
              const double score_cutoff,
              vector<vector<GenomicRegion> > &domains,
              const size_t undef_domain_cutoff);

void
build_domains(const vector<vector<SimpleGenomicRegion> > &bins,
              const vector<vector<size_t> > &classes,
              const vector<vector<double> > &scores, // posterior score of classes[i]
              const double score_cutoff,
              vector<vector<GenomicRegion> > &domains,
              const size_t undef_domain_cutoff);

void
pick_domains(const vector<vector<SimpleGenomicRegion> > &bins,
             const vector<vector<double> > &read_counts,
             const vector<vector<double> > &scales,
             const vector<Distro> &distros,
             vector<vector<GenomicRegion> > &domains,
             const double p_value);

void
pick_domains(const vector<vector<SimpleGenomicRegion> > &bins,
             const vector<vector<double> > &read_counts,
             const vector<vector<double> > &scales,
             const vector<SplitDistro> &distros,
             vector<vector<GenomicRegion> > &domains,
             const double cdf_cutoff);

void
pick_domains(const vector<vector<SimpleGenomicRegion> > &bins,
             const vector<vector<double> > &read_counts,
             const vector<vector<double> > &read_counts_a,
             const vector<vector<double> > &read_counts_b,
             const vector<vector<double> > &scales,
             const vector<SplitDistro> &distros,
             vector<vector<GenomicRegion> > &domains,
             const double cdf_cutoff);

void
pick_domains_3s(const vector<vector<SimpleGenomicRegion> > &bins,
                const vector<vector<double> > &read_counts,
                const vector<vector<double> > &scales,
                const vector<SplitDistro> &distros,
                vector<vector<GenomicRegion> > &domains,
                const double cdf_cutoff);

void
pick_domains_3s(const vector<vector<SimpleGenomicRegion> > &bins,
                const vector<vector<double> > &read_counts,
                const vector<vector<double> > &read_counts_a,
                const vector<vector<double> > &read_counts_b,
                const vector<vector<double> > &scales,
                const vector<SplitDistro> &distros,
                vector<vector<GenomicRegion> > &domains,
                const double cdf_cutoff);


