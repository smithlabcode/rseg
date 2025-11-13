/* Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
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

#include "rseg_utils.hpp"
#include "Distro.hpp"
#include "Interval.hpp"
#include "Interval6.hpp"
#include "SplitDistro.hpp"
#include "log_sum_log.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

template <typename T, typename U>
[[nodiscard]] static auto
overlaps(const T &a, const U &b) {
  return a.chrom == b.chrom &&
         std::max(a.start, b.start) < std::min(a.stop, b.stop);
}

void
pick_training_sample(
  const std::vector<double> &read_bins, const std::vector<double> &read_bins_a,
  const std::vector<double> &read_bins_b, const std::vector<double> &scales,
  const std::vector<std::size_t> &reset_points,
  const std::size_t training_sample_size, std::vector<double> &read_bins_sample,
  std::vector<double> &read_bins_a_sample,
  std::vector<double> &read_bins_b_sample, std::vector<double> &scales_sample,
  std::vector<std::size_t> &reset_points_sample) {
  std::vector<std::size_t> reset_points_control;
  for (std::size_t i = 0; i < reset_points.size(); ++i)
    if (reset_points[i] < training_sample_size) {
      reset_points_sample.push_back(reset_points[i]);
    }
    else {
      reset_points_sample.push_back(training_sample_size);
      break;
    }

  const std::size_t sample_sz = reset_points_sample.back();

  std::copy(read_bins.begin(), read_bins.begin() + sample_sz,
            std::back_inserter(read_bins_sample));
  std::copy(read_bins_a.begin(), read_bins_a.begin() + sample_sz,
            std::back_inserter(read_bins_a_sample));
  std::copy(read_bins_b.begin(), read_bins_b.begin() + sample_sz,
            std::back_inserter(read_bins_b_sample));
  std::copy(scales.begin(), scales.begin() + sample_sz,
            std::back_inserter(scales_sample));
}

void
clear_training_sample(std::vector<double> &read_bins_sample,
                      std::vector<double> &read_bins_a_sample,
                      std::vector<double> &read_bins_b_sample,
                      std::vector<double> &scales_sample,
                      std::vector<std::size_t> &reset_points_sample) {
  read_bins_sample.clear();
  read_bins_sample.shrink_to_fit();

  read_bins_a_sample.clear();
  read_bins_a_sample.shrink_to_fit();

  read_bins_b_sample.clear();
  read_bins_b_sample.shrink_to_fit();

  scales_sample.clear();
  scales_sample.shrink_to_fit();

  reset_points_sample.clear();
  reset_points_sample.shrink_to_fit();
}

void
set_transitions(const std::size_t bin_size, const double fg_size,
                const std::vector<double> &mixing,
                std::vector<double> &start_trans,
                std::vector<std::vector<double>> &trans,
                std::vector<double> &end_trans) {
  start_trans.resize(3, 0);
  start_trans[0] = mixing[0];
  start_trans[1] = mixing[1];
  start_trans[2] = mixing[2];

  end_trans.resize(3, 0);
  end_trans[0] = 1e-10;
  end_trans[1] = 1e-10;
  end_trans[2] = 1e-10;

  trans.resize(3, std::vector<double>(3, 0));

  double fg_bin_n = fg_size / bin_size;
  if (fg_bin_n <= 1) {
    std::cerr
      << "\n[Warning] RSEG may not work as expected. "
      << "The expected differential domain size is smaller than bin size.\n"
      << '\n';
    fg_bin_n = 2;
  }
  trans[0][0] = 1 - 1.0 / fg_bin_n;
  trans[0][1] = 1.0 / fg_bin_n * mixing[1] / (mixing[1] + mixing[2]);
  trans[0][2] = 1.0 / fg_bin_n * mixing[2] / (mixing[1] + mixing[2]);

  // assumming all fg and bg domain are sperated by middle state regions
  double mid_bin_n = mixing[1] * fg_bin_n / (2 * mixing[0]);
  if (mid_bin_n <= 1) {
    std::cerr << "\n[Warning] RSEG may not work as expected. "
              << "The expected basal domain size is smaller than bin size.\n"
              << '\n';
    mid_bin_n = 2;
  }
  trans[1][1] = 1 - 1.0 / mid_bin_n;
  trans[1][0] = 1.0 / mid_bin_n * mixing[0] / (mixing[0] + mixing[2]);
  trans[1][2] = 1.0 / mid_bin_n * mixing[2] / (mixing[0] + mixing[2]);

  double bg_bin_n = fg_bin_n;  // assuming symmetry
  trans[2][2] = 1 - 1.0 / bg_bin_n;
  trans[2][0] = 1.0 / bg_bin_n * mixing[0] / (mixing[0] + mixing[1]);
  trans[2][1] = 1.0 / bg_bin_n * mixing[1] / (mixing[0] + mixing[1]);

  assert(start_trans[0] > 0 && start_trans[1] > 0 && start_trans[2] > 0);
  assert(end_trans[0] > 0 && end_trans[1] > 0 && end_trans[2] > 0);
  assert(trans[0][0] > 0 && trans[0][1] > 0 && trans[0][2] > 0 &&
         trans[1][0] > 0 && trans[1][1] > 0 && trans[1][2] > 0 &&
         trans[2][0] > 0 && trans[2][1] > 0 && trans[2][2] > 0);
}

void
set_transitions(const std::size_t bin_size, const double fg_size,
                const double mixing, std::vector<double> &start_trans,
                std::vector<std::vector<double>> &trans,
                std::vector<double> &end_trans) {
  start_trans.resize(2, 0);
  start_trans[0] = mixing;
  start_trans[1] = 1 - mixing;

  end_trans.resize(2, 0);
  end_trans[0] = 1e-10;
  end_trans[1] = 1e-10;

  trans.resize(2, std::vector<double>(2, 0));
  // foreground
  double fg_bin_n = fg_size / bin_size;
  if (fg_bin_n <= 1) {
    std::cerr << "[Warning] RSEG may not work as expected. "
              << "The expected foreground domain size is smaller than bin size."
              << '\n';
    fg_bin_n = 2;
  }
  trans[0][1] = 1.0 / fg_bin_n;
  trans[0][0] = 1 - trans[0][1];

  double bg_bin_n = fg_bin_n * (1 - mixing) / mixing;
  if (bg_bin_n <= 1) {
    std::cerr
      << "\n[Warning] RSEG may not work as expected."
      << "The expected background domain size is smaller than bin size.\n"
      << '\n';
    bg_bin_n = 2;
  }
  trans[1][0] = 1.0 / bg_bin_n;
  trans[1][1] = 1 - trans[1][0];

  assert(start_trans[0] > 0 && start_trans[1] > 0);
  assert(end_trans[0] > 0 && end_trans[1] > 0);
  assert(trans[0][1] > 0 && trans[0][0] > 0 && trans[1][0] > 0 &&
         trans[1][1] > 0);
}

void
report_final_values(const std::vector<Distro> &distros,
                    const std::vector<std::vector<double>> &trans) {
  std::cout << "FINAL ESTIMATES" << '\n' << "---------------" << '\n';

  std::cout << "Emission distributions" << '\n';
  for (std::size_t i = 0; i < distros.size(); ++i)
    std::cout << "State " << i << ":\t" << distros[i] << '\n';

  std::cout << "Expected sizes" << '\n';
  for (std::size_t i = 0; i < trans.size(); ++i)
    std::cout << "State " << i << ":\t" << 1 / (1 - trans[i][i]) << '\n';

  std::cout << "Transition probabilities" << '\n';
  for (std::size_t i = 0; i < trans.size(); ++i) {
    for (std::size_t j = 0; j < trans[i].size(); ++j)
      std::cout << trans[i][j] << "\t";
    std::cout << '\n';
  }
}

void
report_final_values(const std::vector<SplitDistro> &distros,
                    const std::vector<std::vector<double>> &trans) {
  std::cout << "FINAL ESTIMATES" << '\n' << "---------------" << '\n';

  std::cout << "Emission distributions" << '\n';
  for (std::size_t i = 0; i < distros.size(); ++i)
    std::cout << "State " << i << ":\t" << distros[i] << '\n';

  std::cout << "Expected sizes" << '\n';
  for (std::size_t i = 0; i < trans.size(); ++i)
    std::cout << "State " << i << ":\t" << 1 / (1 - trans[i][i]) << '\n';

  std::cout << "Transition probabilities" << '\n';
  for (std::size_t i = 0; i < trans.size(); ++i) {
    for (std::size_t j = 0; j < trans[i].size(); ++j)
      std::cout << trans[i][j] << "\t";
    std::cout << '\n';
  }
}

void
write_read_counts_by_bin(
  const std::vector<std::vector<Interval>> &bin_boundaries,
  const std::vector<double> &read_bins, const std::vector<double> &scales,
  const std::vector<bool> &classes, const std::string &outfile) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("cannot open file: " + outfile);
  std::size_t k = 0;
  for (std::size_t i = 0; i < bin_boundaries.size(); ++i)
    for (std::size_t j = 0; j < bin_boundaries[i].size(); ++j) {
      std::println(out, "{}\t{}\t{}\t{}", bin_boundaries[i][j], read_bins[k],
                   scales[k], classes[k]);
      ++k;
    }
}

void
write_read_counts_by_bin(
  const std::vector<std::vector<Interval>> &bin_boundaries,
  const std::vector<double> &read_bins, const std::vector<double> &read_bins_a,
  const std::vector<double> &read_bins_b, const std::vector<bool> &classes,
  const std::string &outfile) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("cannot open file: " + outfile);
  std::size_t k = 0;
  for (std::size_t i = 0; i < bin_boundaries.size(); ++i)
    for (std::size_t j = 0; j < bin_boundaries[i].size(); ++j) {
      std::println(out, "{}\t{}\t{}\t{}\t{}", bin_boundaries[i][j],
                   read_bins[k], read_bins_a[k], read_bins_b[k], classes[k]);
      ++k;
    }
}

void
write_read_counts_by_bin(
  const std::vector<std::vector<Interval>> &bin_boundaries,
  const std::vector<double> &read_bins, const std::vector<double> &read_bins_a,
  const std::vector<double> &read_bins_b,
  const std::vector<std::size_t> &classes, const std::string &outfile) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("cannot open file: " + outfile);

  std::size_t k = 0;
  for (std::size_t i = 0; i < bin_boundaries.size(); ++i)
    for (std::size_t j = 0; j < bin_boundaries[i].size(); ++j) {
      std::println(out, "{}\t{}\t{}\t{}\t{}", bin_boundaries[i][j],
                   read_bins[k], read_bins_a[k], read_bins_b[k], classes[k]);
      ++k;
    }
}

std::string
strip_path_and_bed_suffix(const std::string &full_path) {
  std::size_t start = full_path.find_last_of('/');
  if (start == std::string::npos)
    start = 0;
  else
    ++start;
  std::size_t end = full_path.find_last_of('.');
  if (end == std::string::npos)
    end = full_path.length();
  return full_path.substr(start, end - start);
}

void
write_wigfile(const std::vector<std::vector<double>> &scores,
              const std::vector<std::vector<Interval>> &bin_bounds,
              const std::string &wigfile_name) {
  std::ofstream wigout(wigfile_name);
  if (!wigout)
    throw std::runtime_error("cannot open: " + wigfile_name);
  for (std::size_t i = 0; i < bin_bounds.size(); ++i)
    for (std::size_t j = 0; j < bin_bounds[i].size(); ++j)
      std::println(wigout, "{}\t{}", bin_bounds[i][j], scores[i][j]);
}

void
write_wigfile(const std::vector<double> &fg_scores,
              const std::vector<double> &bg_scores,
              const std::vector<std::vector<Interval>> &bin_bounds,
              const std::string &wigfile_name) {
  std::ofstream wigout(wigfile_name);
  if (!wigout)
    throw std::runtime_error("cannot open: " + wigfile_name);
  std::size_t k = 0;
  for (std::size_t i = 0; i < bin_bounds.size(); ++i)
    for (std::size_t j = 0; j < bin_bounds[i].size(); ++j) {
      std::println(wigout, "{}\t{}\t{}", bin_bounds[i][j], fg_scores[k],
                   bg_scores[k]);
      ++k;
    }
}

void
write_bed_file(const std::vector<std::vector<Interval6>> &regions,
               const std::string &outfile) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open: " + outfile);
  for (const auto &region_set : regions)
    for (const auto &region : region_set)
      std::println(out, "{}", region);
}

// for two-state segmentation
void
build_domains(const std::vector<std::vector<Interval>> &bins,
              const std::vector<std::vector<bool>> &classes,
              // 'scores' is posterior score of classes[i]
              const std::vector<std::vector<double>> &scores,
              const double score_cutoff,
              std::vector<std::vector<Interval6>> &domains,
              const std::size_t undef_domain_cutoff =
                std::numeric_limits<std::size_t>::max()) {
  static const std::size_t BG_LABEL = 0;
  static const std::size_t FG_LABEL = 1;
  static const std::size_t UN_LABEL = 2;

  const auto LABEL_NAMES = std::vector<std::string>{
    std::string{"BACKGROUND"},
    std::string{"ENRICHED"},
    std::string{"UNCONFIDENT"},
  };

  std::vector<std::vector<std::size_t>> labels(classes.size());
  std::vector<std::vector<double>> local_scores(scores);

  // STEP I: Relabeling bins
  for (std::size_t i = 0; i < classes.size(); ++i) {
    labels[i].resize(classes[i].size(), UN_LABEL);

    for (std::size_t j = 0; j < classes[i].size(); ++j)
      if (scores[i][j] >= score_cutoff)
        labels[i][j] = static_cast<std::size_t>(classes[i][j]);
  }

  for (std::size_t i = 0; i < bins.size(); ++i) {
    const int lim = static_cast<int>(std::size(bins[i]));

    for (int j = 1; j < lim; ++j)
      if (labels[i][j] == UN_LABEL && labels[i][j - 1] != UN_LABEL &&
          static_cast<std::size_t>(classes[i][j]) == labels[i][j - 1])
        labels[i][j] = labels[i][j - 1];

    for (int j = lim - 2; j >= 0; --j)
      if (labels[i][j] == UN_LABEL && labels[i][j + 1] != UN_LABEL &&
          static_cast<std::size_t>(classes[i][j]) == labels[i][j + 1])
        labels[i][j] = labels[i][j + 1];

    // deal with undefined regions
    int start = 0;
    while (start < lim) {
      while (start < lim && labels[i][start] != UN_LABEL)
        ++start;

      int end = start + 1;
      while (end < lim && labels[i][end] == UN_LABEL)
        ++end;

      if (start >= lim)
        break;  // no undefined bins in this big region

      if (bins[i][end - 1].stop - bins[i][start].start > undef_domain_cutoff) {
        // size of undefined region is big
        start = end;
        continue;
      }

      //// of undefined region is small
      //  determine the likely state of these undefined bins
      const int fg_n = std::accumulate(std::cbegin(classes[i]) + start,
                                       std::cbegin(classes[i]) + end, 0);
      std::size_t label = (fg_n * 2 >= end - start) ? FG_LABEL : BG_LABEL;

      const std::size_t pre_label =
        (start > 0) ? labels[i][start - 1] : UN_LABEL;
      const std::size_t next_label = (end < lim) ? labels[i][end] : UN_LABEL;
      if (pre_label == next_label && pre_label != UN_LABEL)
        label = pre_label;

      for (int j = start; j < end; ++j) {
        labels[i][j] = label;
        if (static_cast<std::size_t>(classes[i][j]) != label)
          local_scores[i][j] = 1 - scores[i][j];
      }

      start = end;
    }
  }

  // STEP II: Build domains
  domains.resize(bins.size(), std::vector<Interval6>());
  for (std::size_t i = 0; i < bins.size(); ++i) {
    domains[i].push_back(Interval6(bins[i].front().chrom, bins[i].front().start,
                                   bins[i].front().stop, "", 0, '+'));

    double current_score = local_scores[i].front();
    for (std::size_t j = 1; j < labels[i].size(); ++j)
      if (labels[i][j] == labels[i][j - 1])
        current_score += local_scores[i][j];
      else {
        domains[i].back().stop = bins[i][j - 1].stop;
        domains[i].back().name = LABEL_NAMES[labels[i][j - 1]];
        domains[i].back().score = current_score;

        domains[i].push_back(Interval6(bins[i][j].chrom, bins[i][j].start,
                                       bins[i][j].stop, "", 0, '+'));
        current_score = local_scores[i][j];
      }
    domains[i].back().stop = bins[i].back().stop;
    domains[i].back().name = LABEL_NAMES[labels[i].back()];
    domains[i].back().score = current_score;
  }
}

// for three-state segmentation
void
build_domains(const std::vector<std::vector<Interval>> &bins,
              const std::vector<std::vector<std::size_t>> &classes,
              // 'scores' is posterior score of classes[i]
              const std::vector<std::vector<double>> &scores,
              const double score_cutoff,
              std::vector<std::vector<Interval6>> &domains,
              const std::size_t undef_domain_cutoff =
                std::numeric_limits<std::size_t>::max()) {
  // static const std::size_t FG_LABEL = 0;
  // static const std::size_t MG_LABEL = 1;
  // static const std::size_t BG_LABEL = 2;
  static const std::size_t UN_LABEL = 3;

  const auto LABEL_NAMES = std::vector<std::string>{
    std::string{"SAMPLE-I-ENRICHED"},
    std::string{"NO-DIFFERENCE"},
    std::string{"SAMPLE-II-ENRICHED"},
    std::string{"UNCONFIDENT"},
  };

  const auto n_classes = std::size(classes);
  std::vector<std::vector<std::size_t>> labels(n_classes);
  std::vector<std::vector<double>> local_scores(scores);

  // STEP I: Relabeling bins
  for (std::size_t i = 0; i < n_classes; ++i) {
    labels[i].resize(classes[i].size(), UN_LABEL);
    for (std::size_t j = 0; j < std::size(classes[i]); ++j)
      if (scores[i][j] >= score_cutoff)
        labels[i][j] = classes[i][j];
  }

  const auto n_bins = std::size(bins);
  for (std::size_t i = 0; i < n_bins; ++i) {
    const int lim = static_cast<int>(std::size(bins[i]));

    for (int j = 1; j < lim; ++j)
      if (labels[i][j] == UN_LABEL && labels[i][j - 1] != UN_LABEL &&
          classes[i][j] == labels[i][j - 1])
        labels[i][j] = classes[i][j];

    for (int j = lim - 2; j >= 0; --j)
      if (labels[i][j] == UN_LABEL && labels[i][j + 1] != UN_LABEL &&
          classes[i][j] == labels[i][j + 1])
        labels[i][j] = classes[i][j];

    // deal with undefined regions
    int start = 0;
    while (start < lim) {
      while (start < lim && labels[i][start] != UN_LABEL)
        ++start;

      int end = start + 1;
      while (end < lim && labels[i][end] == UN_LABEL)
        ++end;

      if (start >= lim)
        break;  // no undefined bins in this big region

      // size of undefined region is big
      if (bins[i][end - 1].stop - bins[i][start].start > undef_domain_cutoff) {
        start = end;
        continue;
      }

      //// of undefined region is small
      //  determine the likely state of these undefined bins
      std::vector<std::size_t> label_bin_nums(3, 0);
      for (int j = start; j < end; ++j)
        ++label_bin_nums[classes[i][j]];

      std::size_t label = std::max_element(std::cbegin(label_bin_nums),
                                           std::cend(label_bin_nums)) -
                          std::cbegin(label_bin_nums);
      const std::size_t pre_label = start > 0 ? labels[i][start - 1] : UN_LABEL;
      const std::size_t next_label = end < lim ? labels[i][end] : UN_LABEL;

      if (pre_label == next_label && pre_label != UN_LABEL)
        label = pre_label;

      for (int j = start; j < end; ++j) {
        labels[i][j] = label;
        if (classes[i][j] != label)
          local_scores[i][j] = (1.0 - scores[i][j]) / 2.0;
      }
      start = end;
    }
  }

  // STEP II: Build domains
  domains.resize(n_bins);
  for (std::size_t i = 0; i < n_bins; ++i) {
    domains[i].push_back(Interval6(bins[i].front().chrom, bins[i].front().start,
                                   bins[i].front().stop, "", 0, '+'));

    double current_score = local_scores[i].front();
    for (std::size_t j = 1; j < labels[i].size(); ++j)
      if (labels[i][j] == labels[i][j - 1])
        current_score += local_scores[i][j];
      else {
        domains[i].back().stop = bins[i][j - 1].stop;
        domains[i].back().name = LABEL_NAMES[labels[i][j - 1]];
        domains[i].back().score = current_score;

        domains[i].push_back(Interval6(bins[i][j].chrom, bins[i][j].start,
                                       bins[i][j].stop, "", 0, '+'));
        current_score = local_scores[i][j];
      }
    domains[i].back().stop = bins[i].back().stop;
    domains[i].back().name = LABEL_NAMES[labels[i].back()];
    domains[i].back().score = current_score;
  }
}

void
pick_domains(const std::vector<std::vector<Interval>> &bins,
             const std::vector<std::vector<double>> &read_counts,
             const std::vector<std::vector<double>> &scales,
             const std::vector<Distro> &distros,
             std::vector<std::vector<Interval6>> &domains,
             const double cdf_cutoff = 0.4) {
  static const std::size_t BG_LABEL = 0;
  static const std::size_t FG_LABEL = 1;
  static const std::size_t UN_LABEL = 2;

  const auto LABEL_NAMES = std::vector<std::string>{
    std::string{"BACKGROUND"},
    std::string{"ENRICHED"},
    std::string{"UNCONFIDENT"},
  };

  std::vector<std::vector<double>> domain_means;

  double max_count(0);

  domain_means.resize(domains.size());
  for (std::size_t i = 0; i < read_counts.size(); ++i) {
    domain_means[i].resize(domains[i].size());

    std::size_t k = 0;
    double domain_read_count = 0;
    double domain_size = 0;
    for (std::size_t j = 0; j < read_counts[i].size(); ++j) {
      if (!overlaps(domains[i][k], bins[i][j])) {
        domain_means[i][k] = domain_read_count / domain_size;
        domain_read_count = 0;
        domain_size = 0;
        ++k;
      }
      domain_read_count += read_counts[i][j];
      domain_size += scales[i][j];

      if (read_counts[i][j] > max_count)
        max_count = read_counts[i][j];
    }
    domain_means[i][k] = domain_read_count / domain_size;
  }

  std::vector<double> fg_cdfs(max_count + 1);
  std::vector<double> bg_cdfs(max_count + 1);

  fg_cdfs[0] = std::exp(distros.front().log_likelihood(0));
  bg_cdfs[0] = std::exp(distros.back().log_likelihood(0));

  for (std::size_t i = 1; i < fg_cdfs.size(); ++i) {
    fg_cdfs[i] = fg_cdfs[i - 1] + std::exp(distros.front().log_likelihood(i));
    bg_cdfs[i] = bg_cdfs[i - 1] + std::exp(distros.back().log_likelihood(i));
  }

  for (std::size_t i = 0; i < domains.size(); ++i)
    for (std::size_t j = 0; j < domains[i].size(); ++j) {
      auto name = domains[i][j].name;
      const std::size_t c = std::floor(domain_means[i][j]);
      if (name == LABEL_NAMES[FG_LABEL] && fg_cdfs[c] < cdf_cutoff)
        name = LABEL_NAMES[UN_LABEL];
      if (name == LABEL_NAMES[BG_LABEL] && 1 - bg_cdfs[c] < cdf_cutoff)
        name = LABEL_NAMES[UN_LABEL];
      domains[i][j].name = name + "\t" + std::to_string(domain_means[i][j]);
    }

  const auto is_enriched = [](const auto &x) {
    return x.name.find("ENRICHED") != std::string::npos;
  };

  for (std::size_t i = 0; i < std::size(domains); ++i)
    domains[i].erase(std::stable_partition(std::begin(domains[i]),
                                           std::end(domains[i]), is_enriched),
                     std::end(domains[i]));
}

void
pick_domains(const std::vector<std::vector<Interval>> &bins,
             const std::vector<std::vector<double>> &read_counts,
             const std::vector<std::vector<double>> &scales,
             const std::vector<SplitDistro> &distros,
             std::vector<std::vector<Interval6>> &domains,
             const double cdf_cutoff = 0.2) {
  static const std::size_t FG_LABEL = 1;
  static const std::size_t BG_LABEL = 0;
  static const std::size_t UN_LABEL = 2;

  std::vector<std::string> LABEL_NAMES(3);
  LABEL_NAMES[FG_LABEL] = "ENRICHED";
  LABEL_NAMES[BG_LABEL] = "BACKGROUND";
  LABEL_NAMES[UN_LABEL] = "UNCONFIDENT";

  std::vector<std::vector<double>> domain_means;

  double max_count{};
  double min_count{std::numeric_limits<double>::max()};

  domain_means.resize(std::size(domains));
  for (std::size_t i = 0; i < std::size(read_counts); ++i) {
    domain_means[i].resize(std::size(domains[i]));

    std::size_t k = 0;
    double domain_read_count = 0;
    double domain_size = 0;
    for (std::size_t j = 0; j < read_counts[i].size(); ++j) {
      if (!overlaps(domains[i][k], bins[i][j])) {
        domain_means[i][k] = domain_read_count / domain_size;
        domain_read_count = 0;
        domain_size = 0;
        ++k;
      }
      domain_read_count += read_counts[i][j];
      domain_size += scales[i][j];

      max_count = std::max(max_count, read_counts[i][j]);
      min_count = std::min(min_count, read_counts[i][j]);
    }
    domain_means[i][k] = domain_read_count / domain_size;
  }

  const double offset = min_count;

  std::vector<double> fg_cdfs(
    static_cast<std::size_t>(max_count - min_count + 1)),
    bg_cdfs(static_cast<std::size_t>(max_count - min_count + 1));
  fg_cdfs[0] = std::exp(distros.front().log_likelihood(0 + offset));
  bg_cdfs[0] = std::exp(distros.back().log_likelihood(0 + offset));
  for (std::size_t i = 1; i < fg_cdfs.size(); ++i) {
    fg_cdfs[i] =
      fg_cdfs[i - 1] + std::exp(distros.front().log_likelihood(i + offset));
    bg_cdfs[i] =
      bg_cdfs[i - 1] + std::exp(distros.back().log_likelihood(i + offset));
  }

  for (std::size_t i = 0; i < domains.size(); ++i)
    for (std::size_t j = 0; j < domains[i].size(); ++j) {
      auto name = domains[i][j].name;
      const std::size_t c = std::floor(domain_means[i][j]);
      if (name == LABEL_NAMES[FG_LABEL] &&
          fg_cdfs[static_cast<std::size_t>(c - offset)] < cdf_cutoff)
        name = LABEL_NAMES[UN_LABEL];
      if (name == LABEL_NAMES[BG_LABEL] &&
          1.0 - bg_cdfs[c - offset] < cdf_cutoff)
        name = LABEL_NAMES[UN_LABEL];
      domains[i][j].name = name + "\t" + std::to_string(domain_means[i][j]);
    }

  const auto is_enriched = [](const auto &x) {
    return x.name.find("ENRICHED") != std::string::npos;
  };

  for (std::size_t i = 0; i < std::size(domains); ++i)
    domains[i].erase(std::stable_partition(std::begin(domains[i]),
                                           std::end(domains[i]), is_enriched),
                     std::end(domains[i]));
}

void
pick_domains_3s(const std::vector<std::vector<Interval>> &bins,
                const std::vector<std::vector<double>> &read_counts,
                const std::vector<std::vector<double>> &scales,
                const std::vector<SplitDistro> &distros,
                std::vector<std::vector<Interval6>> &domains,
                const double cdf_cutoff = 0.4) {
  static const std::size_t FG_LABEL = 0;
  // static const std::size_t MG_LABEL = 1;
  static const std::size_t BG_LABEL = 2;
  static const std::size_t UN_LABEL = 3;

  const auto LABEL_NAMES = std::vector<std::string>{
    std::string{"SAMPLE-I-ENRICHED"},
    std::string{"NO-DIFFERENCE"},
    std::string{"SAMPLE-II-ENRICHED"},
    std::string{"UNCONFIDENT"},
  };

  std::vector<std::vector<double>> domain_means;

  double max_count{};
  double min_count{std::numeric_limits<double>::max()};

  domain_means.resize(domains.size());
  for (std::size_t i = 0; i < std::size(read_counts); ++i) {
    domain_means[i].resize(domains[i].size());

    std::size_t k = 0;
    double domain_read_count = 0;
    double domain_size = 0;
    for (std::size_t j = 0; j < read_counts[i].size(); ++j) {
      if (!overlaps(domains[i][k], bins[i][j])) {
        domain_means[i][k] = domain_read_count / domain_size;
        domain_read_count = 0;
        domain_size = 0;
        ++k;
      }
      domain_read_count += read_counts[i][j];
      domain_size += scales[i][j];

      if (read_counts[i][j] > max_count)
        max_count = read_counts[i][j];

      if (read_counts[i][j] < min_count)
        min_count = read_counts[i][j];
    }
    domain_means[i][k] = domain_read_count / domain_size;
  }

  const double offset = min_count;

  std::vector<double> fg_cdfs(max_count - min_count + 1);
  std::vector<double> bg_cdfs(max_count - min_count + 1);

  fg_cdfs[0] = std::exp(distros.front().log_likelihood(0 + offset));
  bg_cdfs[0] = std::exp(distros.back().log_likelihood(0 + offset));
  for (std::size_t i = 1; i < fg_cdfs.size(); ++i) {
    fg_cdfs[i] =
      fg_cdfs[i - 1] + std::exp(distros.front().log_likelihood(i + offset));
    bg_cdfs[i] =
      bg_cdfs[i - 1] + std::exp(distros.back().log_likelihood(i + offset));
  }

  for (std::size_t i = 0; i < domains.size(); ++i)
    for (std::size_t j = 0; j < domains[i].size(); ++j) {
      auto name = domains[i][j].name;
      const double cf = std::floor(domain_means[i][j]);
      const double cc = std::ceil(domain_means[i][j]);

      if (name == LABEL_NAMES[FG_LABEL] && cc < distros.back().get_mean())
        name = LABEL_NAMES[BG_LABEL];

      if (name == LABEL_NAMES[BG_LABEL] && cf > distros.front().get_mean())
        name = LABEL_NAMES[FG_LABEL];

      if (name == LABEL_NAMES[FG_LABEL] &&
          fg_cdfs[static_cast<std::size_t>(cf - offset)] < cdf_cutoff)
        name = LABEL_NAMES[UN_LABEL];

      if (name == LABEL_NAMES[BG_LABEL] &&
          1 - bg_cdfs[static_cast<std::size_t>(cc - offset)] < cdf_cutoff)
        name = LABEL_NAMES[UN_LABEL];

      domains[i][j].name = name + "\t" + std::to_string(domain_means[i][j]);
    }

  const auto is_enriched = [](const auto &x) {
    return x.name.find("ENRICHED") != std::string::npos;
  };

  for (std::size_t i = 0; i < std::size(domains); ++i)
    domains[i].erase(std::stable_partition(std::begin(domains[i]),
                                           std::end(domains[i]), is_enriched),
                     std::end(domains[i]));
}
