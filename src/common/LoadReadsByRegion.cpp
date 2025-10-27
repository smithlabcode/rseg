/* Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
 *
 * Author: Qiang Song and Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this software; if not, write to the Free Software Foundation, Inc., 51
 * Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "LoadReadsByRegion.hpp"
#include "GenomicRegion.hpp"
#include "SortGenomicRegion.hpp"
#include "bam_record_utils.hpp"

#include <bamxx.hpp>
#include <smithlab_utils.hpp>

#include <cstdint>
#include <fstream>
#include <iomanip>

[[nodiscard]] static GenomicRegion
get_genomic_region(const bamxx::bam_header &hdr, const bamxx::bam_rec &aln) {
  const std::size_t start_pos = get_pos(aln);
  const std::size_t stop_pos = start_pos + get_l_qseq(aln);
  const std::string chrom = sam_hdr_tid2name(hdr, aln);
  const char strand = bam_is_rev(aln) ? '-' : '+';
  return GenomicRegion(chrom, start_pos, stop_pos, "", 0, strand);
}

template <class T>
bool
check_sorted(const std::vector<T> &regions, bool require_unique = false) {
  if (require_unique) {
    for (std::size_t i = 1; i < regions.size(); ++i)
      if (regions[i] <= regions[i - 1])
        return false;
  }
  else
    for (std::size_t i = 1; i < regions.size(); ++i)
      if (regions[i] < regions[i - 1])
        return false;
  return true;
}

/*************************************************
 * This function takes the names of three files (a reads file, a
 * chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBED(const bool VERBOSE, const std::string &chroms_file,
                     const std::string &reads_file,
                     const std::string &deads_file, const std::size_t bin_size,
                     std::vector<SimpleGenomicRegion> &bin_boundaries,
                     std::vector<double> &read_bins,
                     std::vector<double> &nondead_scales,
                     std::vector<std::size_t> &reset_points,
                     const std::size_t FRAGMENT_LEN,
                     const bool REMOVE_JACKPOT) {
  // get the chroms
  if (VERBOSE)
    std::cout << "[LOADING_DATA] chromosomes\n";
  std::vector<SimpleGenomicRegion> chroms;
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions(chroms);

  // Create bins
  reset_points.push_back(0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    const std::string chrom(chroms[i].get_chrom());
    for (std::size_t j = 0; j < chroms[i].get_width(); j += bin_size)
      bin_boundaries.push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    std::cout << "[LOADING_DATA] reads\n";
  read_bins.resize(bin_boundaries.size(), 0);
  std::ifstream in(reads_file);
  std::string line;
  std::int64_t i = 0;
  GenomicRegion prev_gr;
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (REMOVE_JACKPOT && prev_gr.get_start() == gr.get_start() &&
        prev_gr.get_strand() == gr.get_strand() &&
        prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr < prev_gr)
      throw std::runtime_error("reads not sorted in " + reads_file);

    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand()) {
      gr.set_start(gr.get_start() + half_len);
      gr.set_end(gr.get_start() + 1);
    }
    else {
      gr.set_end(gr.get_end() - half_len);
      gr.set_start(gr.get_end() - 1);
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].get_chrom() < gr.get_chrom() ||
            (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
             bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
                      bin_boundaries[i].get_start() >= gr.get_end()))
      --i;

    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins[i];
  }

  // load the dead zones
  if (VERBOSE)
    std::cerr << "[LOADING_DATA] deadzones\n";
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  in.close();
  in.open(deads_file);
  i = 0;
  while (getline(in, line)) {
    SimpleGenomicRegion gr(line);
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           !bin_boundaries[i].overlaps(gr))
      ++i;
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           bin_boundaries[i].overlaps(gr)) {
      const double dead =
        std::min(bin_boundaries[i].get_end(), gr.get_end()) -
        std::max(bin_boundaries[i].get_start(), gr.get_start());
      nondead_scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
}

/*************************************************
 * This function takes the names of four files (reads file a, reads
 * file b, a chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBED(
  const bool VERBOSE, const std::string &chroms_file,
  const std::string &reads_file_a, const std::string &reads_file_b,
  const std::string &deads_file, const std::size_t bin_size,
  std::vector<SimpleGenomicRegion> &bin_boundaries,
  std::vector<double> &read_bins_a, std::vector<double> &read_bins_b,
  std::vector<double> &nondead_scales, std::vector<std::size_t> &reset_points,
  const std::size_t FRAGMENT_LEN, const bool REMOVE_JACKPOT) {
  // get the chroms
  if (VERBOSE)
    std::cout << "[LOADING_DATA] chromosomes\n";
  std::vector<SimpleGenomicRegion> chroms;
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions(chroms);

  // Create bins
  reset_points.push_back(0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    const std::string chrom(chroms[i].get_chrom());
    for (std::size_t j = 0; j < chroms[i].get_width(); j += bin_size)
      bin_boundaries.push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    std::cout << "[LOADING_DATA] reads\n";
  read_bins_a.resize(bin_boundaries.size(), 0);
  std::ifstream in(reads_file_a);
  std::string line;
  std::int64_t i = 0;
  GenomicRegion prev_gr;
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (REMOVE_JACKPOT && prev_gr.get_start() == gr.get_start() &&
        prev_gr.get_strand() == gr.get_strand() &&
        prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr < prev_gr)
      throw std::runtime_error("reads not sorted in " + reads_file_a);
    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand()) {
      gr.set_start(gr.get_start() + half_len);
      gr.set_end(gr.get_start() + 1);
    }
    else {
      gr.set_end(gr.get_end() - half_len);
      gr.set_start(gr.get_end() - 1);
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].get_chrom() < gr.get_chrom() ||
            (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
             bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
                      bin_boundaries[i].get_start() >= gr.get_end()))
      --i;

    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins_a[i];
  }

  // Load the reads, tabulating counts in bins
  read_bins_b.resize(bin_boundaries.size(), 0);
  in.close();
  in.open(reads_file_b);
  i = 0;
  prev_gr = GenomicRegion();
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (REMOVE_JACKPOT && prev_gr.get_start() == gr.get_start() &&
        prev_gr.get_strand() == gr.get_strand() &&
        prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr < prev_gr)
      throw std::runtime_error("reads not sorted in " + reads_file_b);
    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand()) {
      gr.set_start(gr.get_start() + half_len);
      gr.set_end(gr.get_start() + 1);
    }
    else {
      gr.set_end(gr.get_end() - half_len);
      gr.set_start(gr.get_end() - 1);
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].get_chrom() < gr.get_chrom() ||
            (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
             bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
                      bin_boundaries[i].get_start() >= gr.get_end()))
      --i;
    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins_b[i];
  }

  // load the dead zones
  if (VERBOSE)
    std::cout << "[LOADING_DATA] deadzones\n";
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  in.close();
  in.open(deads_file);
  i = 0;
  while (getline(in, line)) {
    SimpleGenomicRegion gr(line);
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           !bin_boundaries[i].overlaps(gr))
      ++i;
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           bin_boundaries[i].overlaps(gr)) {
      const double dead =
        std::min(bin_boundaries[i].get_end(), gr.get_end()) -
        std::max(bin_boundaries[i].get_start(), gr.get_start());
      nondead_scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
}

/*************************************************
 * This function takes the names of three files (a reads file, a
 * chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBAM(const bool VERBOSE, const std::string &chroms_file,
                     const std::string &reads_file,
                     const std::string &deads_file, const std::size_t bin_size,
                     std::vector<SimpleGenomicRegion> &bin_boundaries,
                     std::vector<double> &read_bins,
                     std::vector<double> &nondead_scales,
                     std::vector<std::size_t> &reset_points,
                     const std::size_t FRAGMENT_LEN,
                     const bool REMOVE_JACKPOT) {
  // get the chroms
  if (VERBOSE)
    std::cout << "[LOADING_DATA] chromosomes\n";
  std::vector<SimpleGenomicRegion> chroms;
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions(chroms);

  // Create bins
  reset_points.push_back(0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    const std::string chrom(chroms[i].get_chrom());
    for (std::size_t j = 0; j < chroms[i].get_width(); j += bin_size)
      bin_boundaries.push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    std::cout << "[LOADING_DATA] reads\n";
  read_bins.resize(bin_boundaries.size(), 0);

  bamxx::bam_in hts(reads_file);
  if (!hts)
    throw std::runtime_error("failed to open input file: " + reads_file);
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  std::int64_t i = 0;

  bamxx::bam_rec aln;
  GenomicRegion prev_gr;

  while (hts.read(hdr, aln)) {
    // ADS: skip reads that have no tid -- they are not mapped
    if (get_tid(aln) == -1)
      continue;

    auto gr = get_genomic_region(hdr, aln);

    if (REMOVE_JACKPOT && prev_gr.get_start() == gr.get_start() &&
        prev_gr.get_strand() == gr.get_strand() &&
        prev_gr.get_chrom() == gr.get_chrom())
      continue;

    if (gr.get_chrom() < prev_gr.get_chrom() ||
        gr.get_start() < prev_gr.get_start() ||
        gr.get_end() < prev_gr.get_end() ||
        gr.get_strand() < prev_gr.get_strand())
      throw std::runtime_error("ERROR: reads not sorted in " + reads_file);

    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand()) {
      gr.set_start(gr.get_start() + half_len);
      gr.set_end(gr.get_start() + 1);
    }
    else {
      gr.set_end(gr.get_end() - half_len);
      gr.set_start(gr.get_end() - 1);
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].get_chrom() < gr.get_chrom() ||
            (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
             bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
                      bin_boundaries[i].get_start() >= gr.get_end()))
      --i;
    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins[i];
  }

  // load the dead zones
  if (VERBOSE)
    std::cerr << "[LOADING_DATA] deadzones\n";
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  std::ifstream dead_in(deads_file);
  i = 0;
  std::string line;
  while (getline(dead_in, line)) {
    SimpleGenomicRegion gr(line);
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           !bin_boundaries[i].overlaps(gr))
      ++i;
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           bin_boundaries[i].overlaps(gr)) {
      const double dead =
        std::min(bin_boundaries[i].get_end(), gr.get_end()) -
        std::max(bin_boundaries[i].get_start(), gr.get_start());
      nondead_scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
}

/*************************************************
 * This function takes the names of four files (reads file a, reads
 * file b, a chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBAM(
  const bool VERBOSE, const std::string &chroms_file,
  const std::string &reads_file_a, const std::string &reads_file_b,
  const std::string &deads_file, const std::size_t bin_size,
  std::vector<SimpleGenomicRegion> &bin_boundaries,
  std::vector<double> &read_bins_a, std::vector<double> &read_bins_b,
  std::vector<double> &nondead_scales, std::vector<std::size_t> &reset_points,
  const std::size_t FRAGMENT_LEN, const bool REMOVE_JACKPOT) {
  // get the chroms
  if (VERBOSE)
    std::cout << "[LOADING_DATA] chromosomes\n";
  std::vector<SimpleGenomicRegion> chroms;
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions(chroms);

  // Create bins
  reset_points.push_back(0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    const std::string chrom(chroms[i].get_chrom());
    for (std::size_t j = 0; j < chroms[i].get_width(); j += bin_size)
      bin_boundaries.push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    std::cout << "[LOADING_DATA] reads\n";
  read_bins_a.resize(bin_boundaries.size(), 0);

  bamxx::bam_in hts_a(reads_file_a);
  if (!hts_a)
    throw std::runtime_error("failed to open input file: " + reads_file_a);
  bamxx::bam_header hdr_a(hts_a);
  if (!hdr_a)
    throw std::runtime_error("failed to read header");

  std::int64_t i = 0;
  bamxx::bam_rec aln;

  GenomicRegion prev_gr;

  while (hts_a.read(hdr_a, aln)) {
    // ADS: skip reads that have no tid -- they are not mapped
    if (get_tid(aln) == -1)
      continue;

    auto gr = get_genomic_region(hdr_a, aln);

    if (REMOVE_JACKPOT && prev_gr.get_start() == gr.get_start() &&
        prev_gr.get_strand() == gr.get_strand() &&
        prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr.get_chrom() < prev_gr.get_chrom() ||
        gr.get_start() < prev_gr.get_start() ||
        gr.get_end() < prev_gr.get_end() ||
        gr.get_strand() < prev_gr.get_strand())
      throw std::runtime_error("ERROR: reads not sorted in " + reads_file_a);

    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand()) {
      gr.set_start(gr.get_start() + half_len);
      gr.set_end(gr.get_start() + 1);
    }
    else {
      gr.set_end(gr.get_end() - half_len);
      gr.set_start(gr.get_end() - 1);
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].get_chrom() < gr.get_chrom() ||
            (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
             bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
                      bin_boundaries[i].get_start() >= gr.get_end()))
      --i;
    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins_a[i];
  }

  // Load the reads, tabulating counts in bins
  read_bins_b.resize(bin_boundaries.size(), 0);

  bamxx::bam_in hts_b(reads_file_b);
  if (!hts_b)
    throw std::runtime_error("failed to open input file: " + reads_file_b);
  bamxx::bam_header hdr_b(hts_b);
  if (!hdr_b)
    throw std::runtime_error("failed to read header");

  i = 0;
  prev_gr = GenomicRegion();
  while (hts_b.read(hdr_b, aln)) {
    // ADS: skip reads that have no tid -- they are not mapped
    if (get_tid(aln) == -1)
      continue;

    auto gr = get_genomic_region(hdr_b, aln);

    if (REMOVE_JACKPOT && prev_gr.get_start() == gr.get_start() &&
        prev_gr.get_strand() == gr.get_strand() &&
        prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr.get_chrom() < prev_gr.get_chrom() ||
        gr.get_start() < prev_gr.get_start() ||
        gr.get_end() < prev_gr.get_end() ||
        gr.get_strand() < prev_gr.get_strand())
      throw std::runtime_error("reads not sorted in " + reads_file_b);

    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand()) {
      gr.set_start(gr.get_start() + half_len);
      gr.set_end(gr.get_start() + 1);
    }
    else {
      gr.set_end(gr.get_end() - half_len);
      gr.set_start(gr.get_end() - 1);
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].get_chrom() < gr.get_chrom() ||
            (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
             bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
                      bin_boundaries[i].get_start() >= gr.get_end()))
      --i;
    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins_b[i];
  }

  // load the dead zones
  if (VERBOSE)
    std::cout << "[LOADING_DATA] deadzones\n";
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  std::ifstream dead_in(deads_file);
  i = 0;
  std::string line;
  while (getline(dead_in, line)) {
    SimpleGenomicRegion gr(line);
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           !bin_boundaries[i].overlaps(gr))
      ++i;
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           bin_boundaries[i].overlaps(gr)) {
      const double dead =
        std::min(bin_boundaries[i].get_end(), gr.get_end()) -
        std::max(bin_boundaries[i].get_start(), gr.get_start());
      nondead_scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
}

void
LoadReadsByRegion(const bool VERBOSE, const std::string &chroms_file,
                  const std::string &reads_file, const std::string &deads_file,
                  const std::size_t bin_size,
                  std::vector<SimpleGenomicRegion> &bin_boundaries,
                  std::vector<double> &read_bins,
                  std::vector<double> &nondead_scales,
                  std::vector<std::size_t> &reset_points,
                  const std::size_t FRAGMENT_LEN, const bool BAM_FORMAT,
                  const bool REMOVE_JACKPOT) {
  if (BAM_FORMAT)
    LoadReadsByRegionBAM(VERBOSE, chroms_file, reads_file, deads_file, bin_size,
                         bin_boundaries, read_bins, nondead_scales,
                         reset_points, FRAGMENT_LEN, REMOVE_JACKPOT);
  else
    LoadReadsByRegionBED(VERBOSE, chroms_file, reads_file, deads_file, bin_size,
                         bin_boundaries, read_bins, nondead_scales,
                         reset_points, FRAGMENT_LEN, REMOVE_JACKPOT);
}

void
LoadReadsByRegion(const bool VERBOSE, const std::string &chroms_file,
                  const std::string &reads_file_a,
                  const std::string &reads_file_b,
                  const std::string &deads_file, const std::size_t bin_size,
                  std::vector<SimpleGenomicRegion> &bin_boundaries,
                  std::vector<double> &read_bins_a,
                  std::vector<double> &read_bins_b,
                  std::vector<double> &nondead_scales,
                  std::vector<std::size_t> &reset_points,
                  const std::size_t FRAGMENT_LEN, const bool BAM_FORMAT,
                  const bool REMOVE_JACKPOT) {
  if (BAM_FORMAT)
    LoadReadsByRegionBAM(VERBOSE, chroms_file, reads_file_a, reads_file_b,
                         deads_file, bin_size, bin_boundaries, read_bins_a,
                         read_bins_b, nondead_scales, reset_points,
                         FRAGMENT_LEN, REMOVE_JACKPOT);
  else
    LoadReadsByRegionBED(VERBOSE, chroms_file, reads_file_a, reads_file_b,
                         deads_file, bin_size, bin_boundaries, read_bins_a,
                         read_bins_b, nondead_scales, reset_points,
                         FRAGMENT_LEN, REMOVE_JACKPOT);
}
