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

#include "LoadReadsByRegion.hpp"
#include "SortGenomicRegion.hpp"

#include "smithlab_utils.hpp"

#include <iomanip>
#include <fstream>

using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::setw;
using std::min;
using std::max;

static void
BinDeads(const vector<SimpleGenomicRegion> &deads,
	 const size_t region_start, const size_t region_end,
	 const size_t bin_size, vector<double> &bins) {
  bins.clear();
  size_t dead_idx = 0;
  for (size_t i = region_start; i < region_end; i += bin_size) {
    size_t counts = 0;
    while (dead_idx < deads.size() && deads[dead_idx].get_end() < i + bin_size) {
      counts += deads[dead_idx].get_end() - max(deads[dead_idx].get_start(), i);
      ++dead_idx;
    }
    if (dead_idx < deads.size() && deads[dead_idx].get_start() < i + bin_size)
      counts += ((i + bin_size) - max(deads[dead_idx].get_start(), i));
    assert(counts <= bin_size);
    bins.push_back(counts);
  }
}

static void
BinReadsCorrectDeadZones(const vector<SimpleGenomicRegion> &dead_zones,
			 const size_t start, const size_t end,
			 const size_t bin_size, vector<double> &dead_scales) {
  vector<double> dead_bins;
  BinDeads(dead_zones, start, end, bin_size, dead_bins);
  dead_scales.clear();
  for (size_t i = 0; i < dead_bins.size(); ++i)
    dead_scales.push_back(1 - dead_bins[i]/bin_size);
}

static void
BinReadsCorrectDeadZones(const vector<vector<SimpleGenomicRegion> > &dead_zones,
			 const vector<SimpleGenomicRegion> &regions,
			 const size_t bin_size,
			 vector< vector<double> > &dead_scales) {
  dead_scales.resize(regions.size());
  for (size_t i = 0; i < regions.size(); ++i)
    BinReadsCorrectDeadZones(dead_zones[i], regions[i].get_start(), regions[i].get_end(),
			     bin_size, dead_scales[i]);
}

static void
separate_regions(const vector<SimpleGenomicRegion> &big_regions,
		 const vector<SimpleGenomicRegion> &regions, 
		 vector<vector<SimpleGenomicRegion> > &sep_regions) {
  size_t rr_id = 0;
  const size_t n_regions = regions.size();
  const size_t n_big_regions = big_regions.size();
  sep_regions.resize(n_big_regions);
  for (size_t i = 0; i < n_big_regions; ++i) {
    const string current_chrom(big_regions[i].get_chrom());
    const size_t current_start = big_regions[i].get_start();
    const size_t current_end = big_regions[i].get_end();
    while (rr_id < n_regions &&
	   (regions[rr_id].get_chrom() < current_chrom ||
	    (regions[rr_id].get_chrom() == current_chrom &&
	     regions[rr_id].get_start() < current_start)))
      ++rr_id;
    while (rr_id < n_regions &&
	   (regions[rr_id].get_chrom() == current_chrom &&
	    regions[rr_id].get_start() < current_end)) {
      sep_regions[i].push_back(regions[rr_id]);
      ++rr_id;
    }
  }
}

template <class T>
static void
collapse_read_bins(const vector<vector<T> > &tmp_read_bins, 
                   vector<T> &read_bins,
                   vector<size_t> &reset_points)
{
    reset_points.clear();
    read_bins.clear();
    reset_points.push_back(0);
    for (size_t i = 0; i < tmp_read_bins.size(); ++i) 
    {
        read_bins.insert(read_bins.end(), 
		       tmp_read_bins[i].begin(), 
                         tmp_read_bins[i].end());
        reset_points.push_back(read_bins.size());
    }
}

/*************************************************
 * This function takes the names of three files (a reads file, a
 * chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegion(const bool VERBOSE,
		  const string &chroms_file, const string &reads_file, 
		  const string &deads_file, 
		  const size_t bin_size, 
		  vector<SimpleGenomicRegion> &chroms, 
		  vector<vector<SimpleGenomicRegion> > &deads,
		  vector<vector<double> > &bins,
		  vector<vector<SimpleGenomicRegion> > &boundaries) {
  
  // get the chroms
  if (VERBOSE)
    cerr << "[LOADING_DATA] chromosomes" << endl;
  chroms.clear();
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions_collapse(chroms);

  // get the dead zones
  vector<SimpleGenomicRegion> raw_deads;
  if (!deads_file.empty()) {
    if (VERBOSE)
      cerr << "[LOADING_DATA] dead zones" << endl;
    ReadBEDFile(deads_file, raw_deads);
    if (!check_sorted(deads, true))
      SortGenomicRegion::sort_regions_collapse(raw_deads);
  }
  separate_regions(chroms, raw_deads, deads);
  
  // Create bins
  bins = vector<vector<double> >(chroms.size());
  boundaries = vector<vector<SimpleGenomicRegion> >(chroms.size());
  for (size_t i = 0; i < chroms.size(); ++i) {
    const string chrom(chroms[i].get_chrom());
    for (size_t j = 0; j < chroms[i].get_width(); j += bin_size) {
      boundaries[i].push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
      bins[i].push_back(0);
    }
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
      cerr << "[LOADING_DATA] reads"  << endl;
  std::ifstream in(reads_file.c_str());
  string line;
  size_t i = 0, j = 0;
  GenomicRegion prev_gr;
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
        continue;
    prev_gr = gr;
    if (gr.pos_strand()) gr.set_end(gr.get_start() + 1);
    else gr.set_start(gr.get_end() - 1);
    while (j < boundaries.size()
           && boundaries[j].front().get_chrom() < gr.get_chrom()) {
      ++j;
      i = 0;
    }
    if (j >= boundaries.size() || !gr.same_chrom(boundaries[j].front()))
      continue;
    assert(gr.same_chrom(boundaries[j].front()));
    while (boundaries[j][i] < gr) ++i;
    ++bins[j][i];
  }
}


void
LoadReadsByRegion(const bool VERBOSE,
		  const std::string &regions_file, 
		  const std::string &reads_file, 
		  const std::string &deads_file,
		  const size_t bin_size, 
          std::vector<SimpleGenomicRegion> &boundaries,
          std::vector<double> &read_bins,
		  std::vector<double> &nondead_scales,
          std::vector<size_t> &reset_points)
{
    vector<SimpleGenomicRegion> chroms;
    vector<vector<SimpleGenomicRegion> > deads_folded;
    vector<vector<double> > bins_folded;
    vector<vector<SimpleGenomicRegion> > boundaries_folded;
    vector<vector<double> > nondead_scales_folded;
    LoadReadsByRegion(VERBOSE, regions_file, reads_file, deads_file, 
                      bin_size, chroms, deads_folded, bins_folded,
                      boundaries_folded);
    BinReadsCorrectDeadZones(
        deads_folded, chroms, bin_size, nondead_scales_folded);
    
    collapse_read_bins(boundaries_folded, boundaries, reset_points);
    collapse_read_bins(bins_folded, read_bins, reset_points);
    collapse_read_bins(nondead_scales_folded, nondead_scales, reset_points);

    assert(boundaries.size() == reset_points.back());
    assert(read_bins.size() == reset_points.back());
    assert(nondead_scales.size() == reset_points.back());
}


/*************************************************
 * This function takes the names of four files (reads file a, reads
 * file b, a chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegion(const bool VERBOSE, const string &chroms_file, 
		  const string &reads_file_a, const string &reads_file_b, 
		  const string &deads_file, const size_t bin_size, 
		  vector<SimpleGenomicRegion> &chroms,
		  vector<vector<SimpleGenomicRegion> > &deads,
		  vector<vector<double> > &bins_a, vector<vector<double> > &bins_b,
		  vector<vector<SimpleGenomicRegion> > &boundaries) {
  
  // get the chroms
  if (VERBOSE)
    cerr << "[LOADING_DATA] chromosomes" << endl;
  chroms.clear();
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions_collapse(chroms);

  // get the dead zones
  vector<SimpleGenomicRegion> raw_deads;
  if (!deads_file.empty()) {
    if (VERBOSE)
      cerr << "[LOADING_DATA] dead zones" << endl;
    ReadBEDFile(deads_file, raw_deads);
    if (!check_sorted(deads, true))
      SortGenomicRegion::sort_regions_collapse(raw_deads);
  }
  separate_regions(chroms, raw_deads, deads);


  // Create bins
  bins_a = vector<vector<double> >(chroms.size());
  bins_b = vector<vector<double> >(chroms.size());
  boundaries = vector<vector<SimpleGenomicRegion> >(chroms.size());
  for (size_t i = 0; i < chroms.size(); ++i) {
    const string chrom(chroms[i].get_chrom());
    for (size_t j = 0; j < chroms[i].get_width(); j += bin_size) {
      boundaries[i].push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
      bins_a[i].push_back(0);
      bins_b[i].push_back(0);
    }
  }


  // Load the reads, tabulating counts in bins
  std::ifstream in(reads_file_a.c_str());
  string line;
  size_t i = 0, j = 0;
  GenomicRegion prev_gr;
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
        continue;
    prev_gr = gr;
    if (gr.pos_strand()) gr.set_end(gr.get_start() + 1);
    else gr.set_start(gr.get_end() - 1);
    while (j < boundaries.size()
           &&boundaries[j].front().get_chrom() < gr.get_chrom()) {
      ++j;
      i = 0;
    }
    if (j >= boundaries.size() || !gr.same_chrom(boundaries[j].front()))
      continue;
    assert(gr.same_chrom(boundaries[j].front()));
    while (boundaries[j][i] < gr) ++i;
    ++bins_a[j][i];
  }

  // Load the reads, tabulating counts in bins
  in.close();
  in.open(reads_file_b.c_str());
  i = 0, j = 0;
  prev_gr = GenomicRegion();
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
        continue;
    prev_gr = gr;
    if (gr.pos_strand()) gr.set_end(gr.get_start() + 1);
    else gr.set_start(gr.get_end() - 1);
    while (j < boundaries.size()
           && boundaries[j].front().get_chrom() < gr.get_chrom()) {
      ++j;
      i = 0;
    }
    if (j >= boundaries.size() || !gr.same_chrom(boundaries[j].front()))
      continue;
    assert(gr.same_chrom(boundaries[j].front()));
    while (boundaries[j][i] < gr) ++i;
    ++bins_b[j][i];
  }
}

void
LoadReadsByRegion(const bool VERBOSE,
		  const std::string &regions_file, 
		  const std::string &reads_file_a, 
		  const std::string &reads_file_b, 
		  const std::string &deads_file,
		  const size_t bin_size, 
          std::vector<SimpleGenomicRegion> &boundaries,
          std::vector<double> &read_bins_a,
          std::vector<double> &read_bins_b,
		  std::vector<double> &nondead_scales,
          std::vector<size_t> &reset_points)
{
    vector<SimpleGenomicRegion> chroms;
    vector<vector<SimpleGenomicRegion> > deads_folded;
    vector<vector<double> > bins_a_folded;
    vector<vector<double> > bins_b_folded;
    vector<vector<SimpleGenomicRegion> > boundaries_folded;
    vector<vector<double> > nondead_scales_folded;
    LoadReadsByRegion(VERBOSE, regions_file, reads_file_a, reads_file_b,
                      deads_file, bin_size, chroms, deads_folded,
                      bins_a_folded, bins_b_folded, boundaries_folded);
    BinReadsCorrectDeadZones(
        deads_folded, chroms, bin_size, nondead_scales_folded);
    
    collapse_read_bins(boundaries_folded, boundaries, reset_points);
    collapse_read_bins(bins_a_folded, read_bins_a, reset_points);
    collapse_read_bins(bins_b_folded, read_bins_b, reset_points);
    collapse_read_bins(nondead_scales_folded, nondead_scales, reset_points);

    assert(boundaries.size() == reset_points.back());
    assert(read_bins_a.size() == reset_points.back());
    assert(read_bins_b.size() == reset_points.back());
    assert(nondead_scales.size() == reset_points.back());
}
