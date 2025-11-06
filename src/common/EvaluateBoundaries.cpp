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

#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "EvaluateBoundaries.hpp"
#include "GenomicRegion.hpp"
#include "log_sum_log.hpp"

void
BoundEval::evaluate(
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const std::vector<std::vector<bool>> &classes,
  const std::vector<std::vector<double>> &scores,
  std::vector<std::vector<GenomicRegion>> &boundaries) const {
  // Separate the class values for each bin into domains of contiguous
  // bins having the same class
  const std::size_t n_regions = classes.size();
  std::vector<std::vector<Domain>> domains(n_regions);
  for (std::size_t i = 0; i < n_regions; ++i) {
    std::size_t prev_end = 0;
    for (std::size_t j = 0; j < classes[i].size(); ++j)
      if (j > 0 && classes[i][j] != classes[i][j - 1]) {
        domains[i].push_back(Domain(scores[i], prev_end, j, classes[i][j - 1]));
        prev_end = j;
      }
    domains[i].push_back(
      Domain(scores[i], prev_end, classes[i].size(), classes[i].back()));
  }

  boundaries.resize(n_regions, std::vector<GenomicRegion>());
  for (std::size_t i = 0; i < domains.size(); ++i) {
    const std::string chrom(bin_bounds[i].front().get_chrom());
    boundaries[i].push_back(
      GenomicRegion(chrom, bin_bounds[i].front().get_start(),
                    bin_bounds[i].front().get_start() + 1, "END", 0, '+'));
    std::size_t offset = 0;
    for (std::size_t j = 0; j < domains[i].size() - 1; ++j) {
      offset += domains[i][j].vals.size();
      const double peak_score = domains[i][j + 1].vals.front();
      const std::string peak_name(
        "B:" + std::to_string(static_cast<std::size_t>(peak_score * 1000)));
      boundaries[i].push_back(GenomicRegion(
        chrom, bin_bounds[i][offset].get_start(),
        bin_bounds[i][offset].get_start() + 1, peak_name, peak_score, '+'));
    }
    boundaries[i].push_back(boundaries[i].back());
  }
}

void
BoundEval::evaluate(
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const std::vector<std::vector<std::size_t>> &classes,
  const std::vector<std::vector<double>> &scores,
  std::vector<std::vector<GenomicRegion>> &boundaries) const {

  // Separate the class values for each bin into domains of contiguous
  // bins having the same class
  const std::size_t n_regions = classes.size();
  std::vector<std::vector<Domain>> domains(n_regions);
  for (std::size_t i = 0; i < n_regions; ++i) {
    std::size_t prev_end = 0;
    for (std::size_t j = 0; j < classes[i].size(); ++j)
      if (j > 0 && classes[i][j] != classes[i][j - 1]) {
        domains[i].push_back(Domain(scores[i], prev_end, j, classes[i][j - 1]));
        prev_end = j;
      }
    domains[i].push_back(
      Domain(scores[i], prev_end, classes[i].size(), classes[i].back()));
  }

  boundaries.resize(n_regions, std::vector<GenomicRegion>());
  for (std::size_t i = 0; i < domains.size(); ++i) {
    const std::string chrom(bin_bounds[i].front().get_chrom());
    boundaries[i].push_back(
      GenomicRegion(chrom, bin_bounds[i].front().get_start(),
                    bin_bounds[i].front().get_start() + 1, "END", 0, '+'));
    std::size_t offset = 0;
    for (std::size_t j = 0; j < domains[i].size() - 1; ++j) {
      offset += domains[i][j].vals.size();
      const double peak_score = domains[i][j + 1].vals.front();
      const std::string peak_name(
        "B:" + std::to_string(static_cast<std::size_t>(peak_score * 1000)));
      boundaries[i].push_back(GenomicRegion(
        chrom, bin_bounds[i][offset].get_start(),
        bin_bounds[i][offset].get_start() + 1, peak_name, peak_score, '+'));
    }
    boundaries[i].push_back(boundaries[i].back());
  }
}

std::string
Domain::tostring() const {
  std::ostringstream ss;
  for (const auto &v : vals)
    ss << v << '\n';
  return ss.str();
}

void
BoundEval::evaluate(
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const std::vector<std::vector<bool>> &classes,
  const std::vector<std::vector<double>> &scores,
  std::vector<std::vector<GenomicRegion>> &boundaries,
  std::vector<std::vector<GenomicRegion>> &boundary_peaks,
  std::vector<std::vector<std::size_t>> &boundary_sizes) const {
  const std::size_t n_regions = classes.size();
  std::vector<std::vector<Domain>> domains(n_regions);
  for (std::size_t i = 0; i < n_regions; ++i) {
    std::size_t prev_end = 0;
    for (std::size_t j = 0; j < classes[i].size(); ++j)
      if (j > 0 && classes[i][j] != classes[i][j - 1]) {
        domains[i].push_back(Domain(scores[i], prev_end, j, classes[i][j - 1]));
        prev_end = j;
      }
    domains[i].push_back(
      Domain(scores[i], prev_end, classes[i].size(), classes[i].back()));
  }

  boundaries.resize(n_regions, std::vector<GenomicRegion>());
  boundary_peaks.resize(n_regions, std::vector<GenomicRegion>());
  boundary_sizes.resize(n_regions, std::vector<std::size_t>());
  for (std::size_t i = 0; i < domains.size(); ++i) {
    const std::string chrom(bin_bounds[i].front().get_chrom());
    boundaries[i].push_back(
      GenomicRegion(chrom, bin_bounds[i].front().get_start(),
                    bin_bounds[i].front().get_start() + 1, "END", 0, '+'));
    boundary_peaks[i].push_back(
      GenomicRegion(chrom, bin_bounds[i].front().get_start(),
                    bin_bounds[i].front().get_start() + 1, "END", 0, '+'));
    boundary_sizes[i].push_back(0);

    std::size_t offset = 0;
    for (std::size_t j = 0; j < domains[i].size() - 1; ++j) {
      double area_under_curve = 0;

      std::size_t left_index = 0;
      for (std::size_t k = 0; k < domains[i][j].vals.size(); ++k) {
        const std::size_t index = domains[i][j].vals.size() - 1 - k;
        if (domains[i][j].vals[index] < 0.01) {
          left_index = k;
          break;
        }
        area_under_curve += domains[i][j].vals[index];
      }

      // TODO: make sure the bin actually at the boundary is always counted,
      // even if it has a low score.
      std::size_t right_index = 0;
      for (std::size_t k = 0; k < domains[i][j + 1].vals.size(); ++k) {
        if (domains[i][j + 1].vals[k] < 0.01) {
          right_index = k;
          break;
        }
        area_under_curve += domains[i][j + 1].vals[k];
      }
      // TODO: see above todo
      // area_under_curve = std::max(0.01, area_under_curve);

      offset += domains[i][j].vals.size();
      const std::size_t bound_start =
        bin_bounds[i][offset - left_index].get_start();
      const std::size_t bound_end =
        bin_bounds[i][offset + right_index].get_end();
      const std::size_t bound_bins = left_index + right_index;
      const double peak_score = domains[i][j + 1].vals.front();
      const double denom = bound_bins > 2 ? bound_bins - 1 : 1;
      const double bound_score =
        peak_score - (area_under_curve - peak_score) / denom;

      const std::string bound_name(
        "B:" + std::to_string(static_cast<std::size_t>(bound_score * 1000)));
      const std::string peak_name(
        "B:" + std::to_string(static_cast<std::size_t>(peak_score * 1000)));
      boundaries[i].push_back(GenomicRegion(chrom, bound_start, bound_end,
                                            bound_name, bound_score, '+'));
      boundary_peaks[i].push_back(GenomicRegion(
        chrom, bin_bounds[i][offset].get_start(),
        bin_bounds[i][offset].get_end(), peak_name, peak_score, '+'));
      boundary_sizes[i].push_back(bound_bins);
    }
    boundaries[i].push_back(
      GenomicRegion(chrom, bin_bounds[i].back().get_start(),
                    bin_bounds[i].back().get_end(), "END", 0, '+'));
    boundary_peaks[i].push_back(boundaries[i].back());
    boundary_sizes[i].push_back(0);
  }
}

void
BoundEval::evaluate(
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const std::vector<std::vector<std::size_t>> &classes,
  const std::vector<std::vector<double>> &scores,
  std::vector<std::vector<GenomicRegion>> &boundaries,
  std::vector<std::vector<GenomicRegion>> &boundary_peaks,
  std::vector<std::vector<std::size_t>> &boundary_sizes) const {

  const std::size_t n_regions = classes.size();
  std::vector<std::vector<Domain>> domains(n_regions);
  for (std::size_t i = 0; i < n_regions; ++i) {
    std::size_t prev_end = 0;
    for (std::size_t j = 0; j < classes[i].size(); ++j)
      if (j > 0 && classes[i][j] != classes[i][j - 1]) {
        domains[i].push_back(Domain(scores[i], prev_end, j, classes[i][j - 1]));
        prev_end = j;
      }
    domains[i].push_back(
      Domain(scores[i], prev_end, classes[i].size(), classes[i].back()));
  }

  boundaries.resize(n_regions, std::vector<GenomicRegion>());
  boundary_peaks.resize(n_regions, std::vector<GenomicRegion>());
  boundary_sizes.resize(n_regions, std::vector<std::size_t>());
  for (std::size_t i = 0; i < domains.size(); ++i) {
    const std::string chrom(bin_bounds[i].front().get_chrom());
    boundaries[i].push_back(
      GenomicRegion(chrom, bin_bounds[i].front().get_start(),
                    bin_bounds[i].front().get_end(), "END", 0, '+'));
    boundary_peaks[i].push_back(
      GenomicRegion(chrom, bin_bounds[i].front().get_start(),
                    bin_bounds[i].front().get_end(), "END", 0, '+'));
    boundary_sizes[i].push_back(0);

    std::size_t offset = 0;
    for (std::size_t j = 0; j < domains[i].size() - 1; ++j) {
      double area_under_curve = 0;

      std::size_t left_index = 0;
      for (std::size_t k = 0; k < domains[i][j].vals.size(); ++k) {
        const std::size_t index = domains[i][j].vals.size() - 1 - k;
        if (domains[i][j].vals[index] < 0.01) {
          left_index = k;
          break;
        }
        area_under_curve += domains[i][j].vals[index];
      }

      // TODO: make sure the bin actually at the boundary is always counted,
      // even if it has a low score.
      std::size_t right_index = 0;
      for (std::size_t k = 0; k < domains[i][j + 1].vals.size(); ++k) {
        if (domains[i][j + 1].vals[k] < 0.01) {
          right_index = k;
          break;
        }
        area_under_curve += domains[i][j + 1].vals[k];
      }
      // TODO: see above todo
      // area_under_curve = max(0.01, area_under_curve);

      offset += domains[i][j].vals.size();
      const std::size_t bound_start =
        bin_bounds[i][offset - left_index].get_start();
      const std::size_t bound_end =
        bin_bounds[i][offset + right_index].get_end();
      const std::size_t bound_bins = left_index + right_index;
      const double peak_score = domains[i][j + 1].vals.front();
      const double denom = bound_bins > 2 ? bound_bins - 1 : 1;
      const double bound_score =
        peak_score - (area_under_curve - peak_score) / denom;

      const std::string bound_name(
        "B:" + std::to_string(static_cast<std::size_t>(bound_score * 1000)));
      const std::string peak_name(
        "B:" + std::to_string(static_cast<std::size_t>(peak_score * 1000)));
      boundaries[i].push_back(GenomicRegion(chrom, bound_start, bound_end,
                                            bound_name, bound_score, '+'));
      boundary_peaks[i].push_back(GenomicRegion(
        chrom, bin_bounds[i][offset].get_start(),
        bin_bounds[i][offset].get_end(), peak_name, peak_score, '+'));
      boundary_sizes[i].push_back(bound_bins);
    }
    boundaries[i].push_back(
      GenomicRegion(chrom, bin_bounds[i].back().get_start(),
                    bin_bounds[i].back().get_end(), "END", 0, '+'));
    boundary_peaks[i].push_back(boundaries[i].back());
    boundary_sizes[i].push_back(0);
  }
}

double
bound_trans_prob(const std::vector<double> &ff_probs,
                 const std::vector<double> &fb_probs,
                 const std::vector<double> &bf_probs,
                 const std::vector<double> &bb_probs, const std::size_t start,
                 const std::size_t end) {
  //     static std::size_t c = 0;

  const std::size_t sz = end - start;

  std::vector<double> log_ff_probs(ff_probs.begin() + start,
                                   ff_probs.begin() + end);
  std::vector<double> log_fb_probs(fb_probs.begin() + start,
                                   fb_probs.begin() + end);
  std::vector<double> log_bf_probs(bf_probs.begin() + start,
                                   bf_probs.begin() + end);
  std::vector<double> log_bb_probs(bb_probs.begin() + start,
                                   bb_probs.begin() + end);

  for (std::size_t i = 0; i < sz; ++i) {
    log_ff_probs[i] = log(log_ff_probs[i]);
    log_fb_probs[i] = log(log_fb_probs[i]);
    log_bf_probs[i] = log(log_bf_probs[i]);
    log_bb_probs[i] = log(log_bb_probs[i]);
  }

  double sum = 0;

  // transition occurs in the start
  double prod = log_fb_probs[0];
  for (std::size_t i = 1; i < sz; ++i)
    prod += log_bb_probs[i] - log_sum_log(log_bb_probs[i], log_bf_probs[i]);

  sum = prod;

  // transition occurs afterwards
  prod -= log_fb_probs[0];
  prod += log_ff_probs[0];
  for (std::size_t i = 1; i < sz; ++i) {
    if (i > 1) {
      prod -= log_fb_probs[i - 1] -
              log_sum_log(log_ff_probs[i - 1], log_fb_probs[i - 1]);
      prod += log_ff_probs[i - 1] -
              log_sum_log(log_ff_probs[i - 1], log_fb_probs[i - 1]);
    }
    prod -= log_bb_probs[i] - log_sum_log(log_bb_probs[i], log_bf_probs[i]);
    prod += log_fb_probs[i] - log_sum_log(log_ff_probs[i], log_fb_probs[i]);
    sum = log_sum_log(sum, prod);
  }

  return exp(sum);
}

void
make_boundary(const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
              const std::vector<double> &trans_scores,
              const std::vector<double> &fg_to_fg_trans_score,
              const std::vector<double> &fg_to_bg_trans_score,
              const std::vector<double> &bg_to_fg_trans_score,
              const std::vector<double> &bg_to_bg_trans_score,
              const std::size_t i, const std::size_t offset,
              const std::size_t left_index, const std::size_t right_index,
              GenomicRegion &bound) {
  const std::string bound_chrom =
    bin_bounds[i][left_index - offset].get_chrom();
  const std::size_t bound_start =
    bin_bounds[i][left_index - offset].get_start();
  const std::size_t bound_end =
    bin_bounds[i][right_index - offset - 1].get_end();

  const std::size_t peak_index =
    std::max_element(trans_scores.begin() + left_index,
                     trans_scores.begin() + right_index) -
    trans_scores.begin();
  const double peak_score = trans_scores[peak_index];
  const std::size_t peak_loc = bin_bounds[i][peak_index - offset].get_start();

  const double bound_score =
    std::max(bound_trans_prob(bg_to_bg_trans_score, bg_to_fg_trans_score,
                              fg_to_bg_trans_score, fg_to_fg_trans_score,
                              left_index, right_index),
             bound_trans_prob(fg_to_fg_trans_score, fg_to_bg_trans_score,
                              bg_to_fg_trans_score, bg_to_bg_trans_score,
                              left_index, right_index));
  const std::string bound_name =
    "B:" + std::to_string(right_index - left_index) + ":" +
    std::to_string(peak_loc) + ":" + std::to_string(peak_score);

  bound.set_chrom(bound_chrom);
  bound.set_start(bound_start);
  bound.set_end(bound_end);
  bound.set_name(bound_name);
  bound.set_score(bound_score);
  bound.set_strand('+');
}

void
BoundEval::evaluate(
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const std::vector<std::size_t> &reset_points,
  const std::vector<bool> &classes, const std::vector<double> &trans_scores,
  const std::vector<double> &fg_to_fg_trans_score,
  const std::vector<double> &fg_to_bg_trans_score,
  const std::vector<double> &bg_to_fg_trans_score,
  const std::vector<double> &bg_to_bg_trans_score, const double cutoff,
  const bool Both_Domain_Ends, std::vector<GenomicRegion> &boundaries) const {
  for (std::size_t i = 0; i < reset_points.size() - 1; ++i) {

    const std::size_t offset = reset_points[i];
    const std::size_t start = reset_points[i];
    const std::size_t end = reset_points[i + 1];

    // find domains and domain ends
    std::vector<std::pair<std::size_t, std::size_t>> domain_ends;
    bool prev_class = false;
    if (classes[start]) {
      prev_class = true;
      domain_ends.push_back(std::make_pair(start, 0));
    }
    for (std::size_t j = start + 1; j < end; ++j)
      if (classes[j] != classes[j - 1]) {
        if (prev_class) {
          domain_ends.back().second = j;
          prev_class = false;
        }
        else {
          domain_ends.push_back(std::make_pair(j, 0));
          prev_class = true;
        }
      }
    if (prev_class)
      domain_ends.back().second = end;

    // build and evaluate boundaries
    for (std::size_t j = 0; j < domain_ends.size(); ++j) {
      std::size_t first_left_index = domain_ends[j].first;
      while (first_left_index > start &&
             trans_scores[first_left_index - 1] > cutoff)
        --first_left_index;

      std::size_t first_right_index = domain_ends[j].first;
      while (first_right_index < end &&
             trans_scores[first_right_index] > cutoff)
        ++first_right_index;

      std::size_t second_left_index = domain_ends[j].second;
      while (second_left_index > start &&
             trans_scores[second_left_index - 1] > cutoff)
        --second_left_index;

      std::size_t second_right_index = domain_ends[j].second;
      while (second_right_index < end &&
             trans_scores[second_right_index] > cutoff)
        ++second_right_index;

      if (Both_Domain_Ends && (first_left_index == first_right_index ||
                               second_left_index == second_right_index))
        continue;

      // output first bound
      if (first_left_index < first_right_index) {
        GenomicRegion bound;
        make_boundary(bin_bounds, trans_scores, fg_to_fg_trans_score,
                      fg_to_bg_trans_score, bg_to_fg_trans_score,
                      bg_to_bg_trans_score, i, offset, first_left_index,
                      first_right_index, bound);
        boundaries.push_back(bound);
      }

      // output second bound
      if (second_left_index < second_right_index) {
        GenomicRegion bound;
        make_boundary(bin_bounds, trans_scores, fg_to_fg_trans_score,
                      fg_to_bg_trans_score, bg_to_fg_trans_score,
                      bg_to_bg_trans_score, i, offset, second_left_index,
                      second_right_index, bound);
        boundaries.push_back(bound);
      }
    }
  }
}

void
BoundEval::evaluate(
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const std::vector<std::size_t> &reset_points,
  // const std::vector<bool> &classes,
  const std::vector<double> &trans_scores,
  const std::vector<double> &fg_to_fg_trans_score,
  const std::vector<double> &fg_to_bg_trans_score,
  const std::vector<double> &bg_to_fg_trans_score,
  const std::vector<double> &bg_to_bg_trans_score, const double cutoff,
  std::vector<GenomicRegion> &boundaries) const {
  for (std::size_t i = 0; i < reset_points.size() - 1; ++i) {
    const std::size_t offset = reset_points[i];
    const std::size_t start = reset_points[i];
    const std::size_t end = reset_points[i + 1];

    std::size_t left_index = start;
    std::size_t right_index = left_index;
    while (left_index < end) {
      while (left_index < end && trans_scores[left_index] < cutoff)
        ++left_index;

      right_index = left_index;
      while (right_index < end && trans_scores[right_index] >= cutoff)
        ++right_index;
      if (left_index < end) {
        const std::string bound_chrom =
          bin_bounds[i][left_index - offset].get_chrom();
        const std::size_t bound_start =
          bin_bounds[i][left_index - offset].get_start();
        const std::size_t bound_end =
          bin_bounds[i][right_index - offset - 1].get_end();

        const std::size_t peak_index =
          std::max_element(trans_scores.begin() + left_index,
                           trans_scores.begin() + right_index) -
          trans_scores.begin();
        const double peak_score = trans_scores[peak_index];
        const std::size_t peak_loc =
          bin_bounds[i][peak_index - offset].get_start();

        const double bound_score =
          std::max(bound_trans_prob(bg_to_bg_trans_score, bg_to_fg_trans_score,
                                    fg_to_bg_trans_score, fg_to_fg_trans_score,
                                    left_index, right_index),
                   bound_trans_prob(fg_to_fg_trans_score, fg_to_bg_trans_score,
                                    bg_to_fg_trans_score, bg_to_bg_trans_score,
                                    left_index, right_index));
        const std::string bound_name =
          "B:" + std::to_string(right_index - left_index) + ":" +
          std::to_string(peak_loc) + ":" + std::to_string(peak_score);
        boundaries.push_back(GenomicRegion(bound_chrom, bound_start, bound_end,
                                           bound_name, bound_score, '+'));
      }
      left_index = right_index;
    }
  }
}

double
bound_trans_prob(
  const std::vector<std::vector<std::vector<double>>> &post_trans,
  const std::size_t from_state, const std::size_t to_state,
  const std::size_t start, const std::size_t end) {
  const std::size_t other_state = 3 - from_state - to_state;

  double sum = 0;

  // transition occurs in the start
  double prod = post_trans[from_state][to_state][start];
  for (std::size_t i = start + 1; i < end; ++i)
    prod *=
      post_trans[to_state][to_state][i] /
      (post_trans[to_state][from_state][i] + post_trans[to_state][to_state][i] +
       post_trans[to_state][other_state][i]);

  sum = prod;

  // transition occurs afterwards
  prod /= post_trans[from_state][to_state][start];
  prod *= post_trans[from_state][from_state][start];
  for (std::size_t i = start + 1; i < end; ++i) {
    if (i > start + 1) {
      prod /= post_trans[from_state][to_state][i - 1] /
              (post_trans[from_state][from_state][i - 1] +
               post_trans[from_state][to_state][i - 1] +
               post_trans[from_state][other_state][i - 1]);

      prod *= post_trans[from_state][from_state][i - 1] /
              (post_trans[from_state][from_state][i - 1] +
               post_trans[from_state][to_state][i - 1] +
               post_trans[from_state][other_state][i - 1]);
    }

    prod /=
      post_trans[to_state][to_state][i] /
      (post_trans[to_state][from_state][i] + post_trans[to_state][to_state][i] +
       post_trans[to_state][other_state][i]);

    prod *= post_trans[from_state][to_state][i] /
            (post_trans[from_state][from_state][i] +
             post_trans[from_state][to_state][i] +
             post_trans[from_state][other_state][i]);

    sum += prod;
  }

  return sum;
}

void
BoundEval::evaluate(
  const std::vector<std::vector<SimpleGenomicRegion>> &bin_bounds,
  const std::vector<std::size_t> &reset_points,
  const std::vector<std::size_t> &classes,
  const std::vector<double> &trans_scores,
  const std::vector<std::vector<std::vector<double>>> &post_trans,
  const double cutoff, std::vector<GenomicRegion> &boundaries) const {
  for (std::size_t i = 0; i < reset_points.size() - 1; ++i) {
    const std::size_t offset = reset_points[i];
    const std::size_t start = reset_points[i];
    const std::size_t end = reset_points[i + 1];

    // find domains and domain ends
    std::vector<std::size_t> domain_ends;
    domain_ends.push_back(start);
    for (std::size_t j = start + 1; j < end; ++j)
      if (classes[j] != classes[j - 1])
        domain_ends.push_back(j);

    domain_ends.push_back(end);

    for (std::size_t j = 0; j < domain_ends.size(); ++j) {
      std::size_t left_index = domain_ends[j];
      while (left_index > start && trans_scores[left_index - 1] > cutoff)
        --left_index;

      std::size_t right_index = domain_ends[j];
      while (right_index < end && trans_scores[right_index] > cutoff)
        ++right_index;

      if (left_index < right_index) {
        const std::string bound_chrom =
          bin_bounds[i][left_index - offset].get_chrom();
        const std::size_t bound_start =
          bin_bounds[i][left_index - offset].get_start();
        const std::size_t bound_end =
          bin_bounds[i][right_index - offset - 1].get_end();

        const std::size_t peak_index =
          std::max_element(trans_scores.begin() + left_index,
                           trans_scores.begin() + right_index) -
          trans_scores.begin();
        const double peak_score = trans_scores[peak_index];
        const std::size_t peak_loc =
          bin_bounds[i][peak_index - offset].get_start();

        const std::size_t state_a = classes[left_index];
        std::size_t state_b, state_c;
        if (state_a == 0) {
          state_b = 1;
          state_c = 2;
        }
        else if (state_a == 1) {
          state_b = 0;
          state_c = 2;
        }
        else  // if (state_a == 2)
        {
          state_b = 0;
          state_c = 1;
        }

        const double bound_score =
          std::max(std::max(bound_trans_prob(post_trans, state_a, state_b,
                                             left_index, right_index),
                            bound_trans_prob(post_trans, state_b, state_a,
                                             left_index, right_index)),
                   std::max(bound_trans_prob(post_trans, state_a, state_c,
                                             left_index, right_index),
                            bound_trans_prob(post_trans, state_c, state_a,
                                             left_index, right_index)));
        const std::string bound_name =
          "B:" + std::to_string(right_index - left_index) + ":" +
          std::to_string(peak_loc) + ":" + std::to_string(peak_score);
        boundaries.push_back(GenomicRegion(bound_chrom, bound_start, bound_end,
                                           bound_name, bound_score, '+'));
      }
    }
  }
}
