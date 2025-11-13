/* Copyright (C) 2025 Andrew D Smith
 * Author: Andrew D. Smith
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

#ifndef EVALUATE_BOUNDARIES_HPP_
#define EVALUATE_BOUNDARIES_HPP_

#include <cstddef>
#include <iterator>
#include <string>
#include <vector>

struct Interval6;
struct Interval;

struct Domain {
  std::vector<double> vals;
  std::size_t state;
  Domain(const std::vector<double> &v, std::size_t start, std::size_t end,
         const bool state) :
    vals{std::cbegin(v) + start, std::cbegin(v) + end},
    state{state} {}
  std::string
  tostring() const;
};

struct BoundEval {
  std::size_t boundary_size;
  double bandwidth;
  void
  evaluate(const std::vector<std::vector<Interval>> &bin_bounds,
           const std::vector<std::vector<bool>> &classes,
           const std::vector<std::vector<double>> &scores,
           std::vector<std::vector<Interval6>> &boundary_scores,
           std::vector<std::vector<Interval6>> &boundary_peaks,
           std::vector<std::vector<std::size_t>> &boundary_sizes) const;
  void
  evaluate(const std::vector<std::vector<Interval>> &bin_bounds,
           const std::vector<std::vector<std::size_t>> &classes,
           const std::vector<std::vector<double>> &scores,
           std::vector<std::vector<Interval6>> &boundary_scores,
           std::vector<std::vector<Interval6>> &boundary_peaks,
           std::vector<std::vector<std::size_t>> &boundary_sizes) const;
  void
  evaluate(const std::vector<std::vector<Interval>> &bin_bounds,
           const std::vector<std::vector<bool>> &classes,
           const std::vector<std::vector<double>> &scores,
           std::vector<std::vector<Interval6>> &boundaries) const;

  void
  evaluate(const std::vector<std::vector<Interval>> &bin_bounds,
           const std::vector<std::vector<std::size_t>> &classes,
           const std::vector<std::vector<double>> &scores,
           std::vector<std::vector<Interval6>> &boundaries) const;

  void
  evaluate(const std::vector<std::vector<Interval>> &bin_bounds,
           const std::vector<std::vector<bool>> &classes,
           const std::vector<std::vector<double>> &scores,
           const std::vector<std::vector<double>> &fg_to_fg_trans_score,
           const std::vector<std::vector<double>> &fg_to_bg_trans_score,
           const std::vector<std::vector<double>> &bg_to_fg_trans_score,
           const std::vector<std::vector<double>> &bg_to_bg_trans_score,
           std::vector<std::vector<Interval6>> &boundaries,
           std::vector<std::vector<Interval6>> &boundary_peaks,
           std::vector<std::vector<std::size_t>> &boundary_sizes) const;

  void
  evaluate(const std::vector<std::vector<Interval>> &bin_bounds,
           const std::vector<std::size_t> &reset_points,
           // const std::vector<bool> &classes,
           const std::vector<double> &trans_scores,
           const std::vector<double> &fg_to_fg_trans_score,
           const std::vector<double> &fg_to_bg_trans_score,
           const std::vector<double> &bg_to_fg_trans_score,
           const std::vector<double> &bg_to_bg_trans_score, const double cutoff,
           std::vector<Interval6> &boundaries) const;

  void
  evaluate(const std::vector<std::vector<Interval>> &bin_bounds,
           const std::vector<std::size_t> &reset_points,
           const std::vector<std::size_t> &classes,
           const std::vector<double> &trans_scores,
           const std::vector<std::vector<std::vector<double>>> &post_trans,
           const double cutoff, std::vector<Interval6> &boundaries) const;
  void
  evaluate(const std::vector<std::vector<Interval>> &bin_bounds,
           const std::vector<std::size_t> &reset_points,
           const std::vector<bool> &classes,
           const std::vector<double> &trans_scores,
           const std::vector<double> &fg_to_fg_trans_score,
           const std::vector<double> &fg_to_bg_trans_score,
           const std::vector<double> &bg_to_fg_trans_score,
           const std::vector<double> &bg_to_bg_trans_score, const double cutoff,
           const bool both_domain_ends,
           std::vector<Interval6> &boundaries) const;
};

#endif  // EVALUATE_BOUNDARIES_HPP_
