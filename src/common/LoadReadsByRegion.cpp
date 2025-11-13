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
#include "Interval.hpp"
#include "Interval6.hpp"
#include "bam_record_utils.hpp"

#include <bamxx.hpp>
// #include <smithlab_utils.hpp>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <print>
#include <ranges>

template <typename T, typename U>
[[nodiscard]] static auto
overlaps(const T &a, const U &b) {
  return a.chrom == b.chrom &&
         std::max(a.start, b.start) < std::min(a.stop, b.stop);
}

template <typename T, typename U>
[[nodiscard]] static auto
contains(const T &bigger, const U &smaller) {
  return bigger.chrom == smaller.chrom && bigger.start <= smaller.start &&
         smaller.stop <= bigger.stop;
}

[[nodiscard]] static Interval6
get_genomic_region(const bamxx::bam_header &hdr, const bamxx::bam_rec &aln) {
  const std::size_t start_pos = get_pos(aln);
  const std::size_t stop_pos = start_pos + get_l_qseq(aln);
  const std::string chrom = sam_hdr_tid2name(hdr, aln);
  const char strand = bam_is_rev(aln) ? '-' : '+';
  return Interval6(chrom, start_pos, stop_pos, "", 0, strand);
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

static auto
get_nondead_scales(const std::size_t bin_size,
                   const std::vector<Interval> &bins,
                   const std::string &deads_file) -> std::vector<double> {
  const std::uint32_t n_bins = std::size(bins);
  std::vector<double> scales(n_bins, 1.0);
  if (deads_file.empty())
    return scales;

  std::ifstream in(deads_file);
  if (!in)
    throw std::runtime_error("failed to open deadzones file: " + deads_file);

  std::uint32_t i = 0;
  std::string line;
  while (getline(in, line)) {
    Interval dead_zone(line);
    while (i < n_bins && !overlaps(bins[i], dead_zone))
      ++i;
    while (i < n_bins && overlaps(bins[i], dead_zone)) {
      const double dead = std::min(bins[i].stop, dead_zone.stop) -
                          std::max(bins[i].start, dead_zone.start);
      scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
  return scales;
}

/*************************************************
 * This function takes the names of three files (a reads file, a
 * chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBED(const bool VERBOSE, const std::string &chroms_file,
                     const std::string &reads_file,
                     const std::string &deads_file, const std::size_t bin_size,
                     std::vector<Interval> &bin_boundaries,
                     std::vector<double> &read_bins,
                     std::vector<double> &nondead_scales,
                     std::vector<std::size_t> &reset_points,
                     const std::size_t FRAGMENT_LEN,
                     const bool REMOVE_JACKPOT) {
  // get the chroms
  if (VERBOSE)
    std::println("[LOADING_DATA] chromosomes");
  auto chroms = read_intervals(chroms_file);
  std::ranges::sort(chroms);

  // Create bins
  reset_points.push_back(0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    const std::string chrom(chroms[i].chrom);
    for (std::size_t j = 0; j < size(chroms[i]); j += bin_size)
      bin_boundaries.push_back(Interval(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    std::println("[LOADING_DATA] reads");
  read_bins.resize(bin_boundaries.size(), 0);
  std::ifstream in(reads_file);
  std::string line;
  std::int64_t i = 0;
  Interval6 prev_gr;
  while (getline(in, line)) {
    Interval6 gr(line);
    if (REMOVE_JACKPOT && prev_gr.start == gr.start &&
        prev_gr.strand == gr.strand && prev_gr.chrom == gr.chrom)
      continue;
    if (gr < prev_gr)
      throw std::runtime_error("reads not sorted in " + reads_file);

    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.strand == '+') {
      gr.start = gr.start + half_len;
      gr.stop = gr.start + 1;
    }
    else {
      gr.stop = gr.stop - half_len;
      gr.start = gr.stop - 1;
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].chrom < gr.chrom ||
            (bin_boundaries[i].chrom == gr.chrom &&
             bin_boundaries[i].stop <= gr.start)))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].chrom == gr.chrom &&
                      bin_boundaries[i].start >= gr.stop))
      --i;

    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !contains(bin_boundaries[i], gr))
      continue;
    ++read_bins[i];
  }

  if (!deads_file.empty())
    std::println(std::cerr, "[LOADING_DATA] deadzones");

  nondead_scales = get_nondead_scales(bin_size, bin_boundaries, deads_file);
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
  std::vector<Interval> &bin_boundaries, std::vector<double> &read_bins_a,
  std::vector<double> &read_bins_b, std::vector<double> &nondead_scales,
  std::vector<std::size_t> &reset_points, const std::size_t FRAGMENT_LEN,
  const bool REMOVE_JACKPOT) {
  // get the chroms
  if (VERBOSE)
    std::println("[LOADING_DATA] chromosomes");
  auto chroms = read_intervals(chroms_file);
  std::ranges::sort(chroms);

  // Create bins
  reset_points.push_back(0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    const std::string chrom(chroms[i].chrom);
    for (std::size_t j = 0; j < size(chroms[i]); j += bin_size)
      bin_boundaries.push_back(Interval(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    std::println("[LOADING_DATA] reads");
  read_bins_a.resize(bin_boundaries.size(), 0);
  std::ifstream in(reads_file_a);
  std::string line;
  std::int64_t i = 0;
  Interval6 prev_gr;
  while (getline(in, line)) {
    Interval6 gr(line);
    if (REMOVE_JACKPOT && prev_gr.start == gr.start &&
        prev_gr.strand == gr.strand && prev_gr.chrom == gr.chrom)
      continue;
    if (gr < prev_gr)
      throw std::runtime_error("reads not sorted in " + reads_file_a);
    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.strand == '+') {
      gr.start = gr.start + half_len;
      gr.stop = gr.start + 1;
    }
    else {
      gr.stop = gr.stop - half_len;
      gr.start = gr.stop - 1;
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].chrom < gr.chrom ||
            (bin_boundaries[i].chrom == gr.chrom &&
             bin_boundaries[i].stop <= gr.start)))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].chrom == gr.chrom &&
                      bin_boundaries[i].start >= gr.stop))
      --i;

    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !contains(bin_boundaries[i], gr))
      continue;
    ++read_bins_a[i];
  }

  // Load the reads, tabulating counts in bins
  read_bins_b.resize(bin_boundaries.size(), 0);
  in.close();
  in.open(reads_file_b);
  i = 0;
  prev_gr = Interval6();
  while (getline(in, line)) {
    Interval6 gr(line);
    if (REMOVE_JACKPOT && prev_gr.start == gr.start &&
        prev_gr.strand == gr.strand && prev_gr.chrom == gr.chrom)
      continue;
    if (gr < prev_gr)
      throw std::runtime_error("reads not sorted in " + reads_file_b);
    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.strand == '+') {
      gr.start = gr.start + half_len;
      gr.stop = gr.start + 1;
    }
    else {
      gr.stop = gr.stop - half_len;
      gr.start = gr.stop - 1;
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].chrom < gr.chrom ||
            (bin_boundaries[i].chrom == gr.chrom &&
             bin_boundaries[i].stop <= gr.start)))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].chrom == gr.chrom &&
                      bin_boundaries[i].start >= gr.stop))
      --i;
    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !contains(bin_boundaries[i], gr))
      continue;
    ++read_bins_b[i];
  }

  if (!deads_file.empty())
    std::println(std::cerr, "[LOADING_DATA] deadzones");

  nondead_scales = get_nondead_scales(bin_size, bin_boundaries, deads_file);
}

/*************************************************
 * This function takes the names of three files (a reads file, a
 * chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBAM(const bool VERBOSE, const std::string &chroms_file,
                     const std::string &reads_file,
                     const std::string &deads_file, const std::size_t bin_size,
                     std::vector<Interval> &bin_boundaries,
                     std::vector<double> &read_bins,
                     std::vector<double> &nondead_scales,
                     std::vector<std::size_t> &reset_points,
                     const std::size_t FRAGMENT_LEN,
                     const bool REMOVE_JACKPOT) {
  // get the chroms
  if (VERBOSE)
    std::println("[LOADING_DATA] chromosomes");
  auto chroms = read_intervals(chroms_file);
  std::ranges::sort(chroms);

  // Create bins
  reset_points.push_back(0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    const std::string chrom(chroms[i].chrom);
    for (std::size_t j = 0; j < size(chroms[i]); j += bin_size)
      bin_boundaries.push_back(Interval(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    std::println("[LOADING_DATA] reads");
  read_bins.resize(bin_boundaries.size(), 0);

  bamxx::bam_in hts(reads_file);
  if (!hts)
    throw std::runtime_error("failed to open input file: " + reads_file);
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  std::int64_t i = 0;

  bamxx::bam_rec aln;
  Interval6 prev_gr;

  while (hts.read(hdr, aln)) {
    // ADS: skip reads that have no tid -- they are not mapped
    if (get_tid(aln) == -1)
      continue;

    auto gr = get_genomic_region(hdr, aln);

    if (REMOVE_JACKPOT && prev_gr.start == gr.start &&
        prev_gr.strand == gr.strand && prev_gr.chrom == gr.chrom)
      continue;

    if (gr.chrom < prev_gr.chrom || gr.start < prev_gr.start ||
        gr.stop < prev_gr.stop || gr.strand < prev_gr.strand)
      throw std::runtime_error("ERROR: reads not sorted in " + reads_file);

    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.strand == '+') {
      gr.start = gr.start + half_len;
      gr.stop = gr.start + 1;
    }
    else {
      gr.stop = gr.stop - half_len;
      gr.start = gr.stop - 1;
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].chrom < gr.chrom ||
            (bin_boundaries[i].chrom == gr.chrom &&
             bin_boundaries[i].stop <= gr.start)))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].chrom == gr.chrom &&
                      bin_boundaries[i].start >= gr.stop))
      --i;
    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !contains(bin_boundaries[i], gr))
      continue;
    ++read_bins[i];
  }

  // load the dead zones
  if (VERBOSE)
    std::println(std::cerr, "[LOADING_DATA] deadzones");
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  std::ifstream dead_in(deads_file);
  i = 0;
  std::string line;
  while (getline(dead_in, line)) {
    Interval gr(line);
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           !overlaps(bin_boundaries[i], gr))
      ++i;
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           overlaps(bin_boundaries[i], gr)) {
      const double dead = std::min(bin_boundaries[i].stop, gr.stop) -
                          std::max(bin_boundaries[i].start, gr.start);
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
  std::vector<Interval> &bin_boundaries, std::vector<double> &read_bins_a,
  std::vector<double> &read_bins_b, std::vector<double> &nondead_scales,
  std::vector<std::size_t> &reset_points, const std::size_t FRAGMENT_LEN,
  const bool REMOVE_JACKPOT) {
  // get the chroms
  if (VERBOSE)
    std::println("[LOADING_DATA] chromosomes");
  auto chroms = read_intervals(chroms_file);
  std::ranges::sort(chroms);

  // Create bins
  reset_points.push_back(0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    const std::string chrom(chroms[i].chrom);
    for (std::size_t j = 0; j < size(chroms[i]); j += bin_size)
      bin_boundaries.push_back(Interval(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    std::println("[LOADING_DATA] reads");
  read_bins_a.resize(bin_boundaries.size(), 0);

  bamxx::bam_in hts_a(reads_file_a);
  if (!hts_a)
    throw std::runtime_error("failed to open input file: " + reads_file_a);
  bamxx::bam_header hdr_a(hts_a);
  if (!hdr_a)
    throw std::runtime_error("failed to read header");

  std::int64_t i = 0;
  bamxx::bam_rec aln;

  Interval6 prev_gr;

  while (hts_a.read(hdr_a, aln)) {
    // ADS: skip reads that have no tid -- they are not mapped
    if (get_tid(aln) == -1)
      continue;

    auto gr = get_genomic_region(hdr_a, aln);

    if (REMOVE_JACKPOT && prev_gr.start == gr.start &&
        prev_gr.strand == gr.strand && prev_gr.chrom == gr.chrom)
      continue;
    if (gr.chrom < prev_gr.chrom || gr.start < prev_gr.start ||
        gr.stop < prev_gr.stop || gr.strand < prev_gr.strand)
      throw std::runtime_error("ERROR: reads not sorted in " + reads_file_a);

    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.strand == '+') {
      gr.start = gr.start + half_len;
      gr.stop = gr.start + 1;
    }
    else {
      gr.stop = gr.stop - half_len;
      gr.start = gr.stop - 1;
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].chrom < gr.chrom ||
            (bin_boundaries[i].chrom == gr.chrom &&
             bin_boundaries[i].stop <= gr.start)))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].chrom == gr.chrom &&
                      bin_boundaries[i].start >= gr.stop))
      --i;
    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !contains(bin_boundaries[i], gr))
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
  prev_gr = Interval6();
  while (hts_b.read(hdr_b, aln)) {
    // ADS: skip reads that have no tid -- they are not mapped
    if (get_tid(aln) == -1)
      continue;

    auto gr = get_genomic_region(hdr_b, aln);

    if (REMOVE_JACKPOT && prev_gr.start == gr.start &&
        prev_gr.strand == gr.strand && prev_gr.chrom == gr.chrom)
      continue;
    if (gr.chrom < prev_gr.chrom || gr.start < prev_gr.start ||
        gr.stop < prev_gr.stop || gr.strand < prev_gr.strand)
      throw std::runtime_error("reads not sorted in " + reads_file_b);

    prev_gr = gr;

    const std::size_t half_len = FRAGMENT_LEN / 2;
    if (gr.strand == '+') {
      gr.start = gr.start + half_len;
      gr.stop = gr.start + 1;
    }
    else {
      gr.stop = gr.stop - half_len;
      gr.start = gr.stop - 1;
    }

    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           (bin_boundaries[i].chrom < gr.chrom ||
            (bin_boundaries[i].chrom == gr.chrom &&
             bin_boundaries[i].stop <= gr.start)))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 && (bin_boundaries[i].chrom == gr.chrom &&
                      bin_boundaries[i].start >= gr.stop))
      --i;
    if (i >= static_cast<std::int64_t>(std::size(bin_boundaries)) ||
        !contains(bin_boundaries[i], gr))
      continue;
    ++read_bins_b[i];
  }

  // load the dead zones
  if (VERBOSE)
    std::println("[LOADING_DATA] deadzones");
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  std::ifstream dead_in(deads_file);
  i = 0;
  std::string line;
  while (getline(dead_in, line)) {
    Interval gr(line);
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           !overlaps(bin_boundaries[i], gr))
      ++i;
    while (i < static_cast<std::int64_t>(std::size(bin_boundaries)) &&
           overlaps(bin_boundaries[i], gr)) {
      const double dead = std::min(bin_boundaries[i].stop, gr.stop) -
                          std::max(bin_boundaries[i].start, gr.start);
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
                  std::vector<Interval> &bin_boundaries,
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
LoadReadsByRegion(
  const bool VERBOSE, const std::string &chroms_file,
  const std::string &reads_file_a, const std::string &reads_file_b,
  const std::string &deads_file, const std::size_t bin_size,
  std::vector<Interval> &bin_boundaries, std::vector<double> &read_bins_a,
  std::vector<double> &read_bins_b, std::vector<double> &nondead_scales,
  std::vector<std::size_t> &reset_points, const std::size_t FRAGMENT_LEN,
  const bool BAM_FORMAT, const bool REMOVE_JACKPOT) {
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
