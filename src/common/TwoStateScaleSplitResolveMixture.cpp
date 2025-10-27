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

#include "TwoStateScaleSplitResolveMixture.hpp"

#include "SplitDistro.hpp"

#include "smithlab_utils.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>

static void
initialize_values(const std::vector<double> &values,
                  const std::vector<double> &values_a,
                  const std::vector<double> &values_b,
                  const std::vector<double> &scales, SplitDistro &fg_distro,
                  SplitDistro &bg_distro) {
  std::vector<double> fg_part_a;
  std::vector<double> fg_part_b;
  std::vector<double> bg_part_a;
  std::vector<double> bg_part_b;
  std::vector<double> fg_scales;
  std::vector<double> bg_scales;

  double val_weight = 0;
  for (std::size_t i = 0; i < values.size(); ++i)
    val_weight += fabs(values[i]);

  assert(val_weight > 0);
  double val_per_part = std::ceil(val_weight / 2);

  std::vector<std::pair<double, std::size_t>> helper;
  for (std::size_t i = 0; i < values.size(); ++i)
    helper.push_back(std::make_pair(values[i], i));
  std::sort(std::begin(helper), std::end(helper));

  double sum = 0;
  std::size_t i = 0;
  for (; i < helper.size() && sum < val_per_part; ++i) {
    fg_part_a.push_back(values_a[helper[i].second]);
    fg_part_b.push_back(values_b[helper[i].second]);
    fg_scales.push_back(scales[helper[i].second]);
    sum += fabs(helper[i].first);
  }
  for (; i < helper.size(); ++i) {
    bg_part_a.push_back(values_a[helper[i].second]);
    bg_part_b.push_back(values_b[helper[i].second]);
    bg_scales.push_back(scales[helper[i].second]);
  }

  // Obtain initial estimates of distribution values
  fg_distro.estimate_params_ml(fg_part_a, fg_part_b, fg_scales);
  bg_distro.estimate_params_ml(bg_part_a, bg_part_b, bg_scales);
}

static double
expectation_step(const std::vector<double> &values,
                 const std::vector<double> &scales, const double mixing,
                 const SplitDistro &fg_distro, const SplitDistro &bg_distro,
                 std::vector<double> &fg_probs, std::vector<double> &bg_probs) {
  double score = 0;

  const double fg_log_mixing = log(mixing);
  assert(std::isfinite(fg_log_mixing));
  const double bg_log_mixing = log(1 - mixing);
  assert(std::isfinite(bg_log_mixing));

  for (std::size_t i = 0; i < std::size(values); ++i) {
    const double fg_part =
      fg_log_mixing + fg_distro.log_likelihood(values[i], scales[i]);
    assert(std::isfinite(fg_part));

    const double bg_part =
      bg_log_mixing + bg_distro.log_likelihood(values[i], scales[i]);
    assert(std::isfinite(fg_part));

    const double denom = fg_part > bg_part
                           ? fg_part + log(1.0 + exp(bg_part - fg_part))
                           : bg_part + log(1.0 + exp(fg_part - bg_part));
    assert(std::isfinite(denom));

    fg_probs[i] = exp(fg_part - denom);
    bg_probs[i] = exp(bg_part - denom);

    score += denom;
  }
  return score;
}

static void
maximization_step(  // const std::vector<double> &values,
                    // clang-format off
   const std::vector<double> &vals_a,
   const std::vector<double> &vals_b,
   const std::vector<double> &scales,
   const std::vector<double> &fg_probs,
   const std::vector<double> &bg_probs,
   double &mixing,
   SplitDistro &fg_distro,
   SplitDistro &bg_distro
                    // clang-format on
) {
  fg_distro.estimate_params_ml(vals_a, vals_b, scales, fg_probs);
  bg_distro.estimate_params_ml(vals_a, vals_b, scales, bg_probs);

  std::vector<double> log_fg_probs(fg_probs);
  std::vector<double> log_bg_probs(bg_probs);
  for (std::size_t i = 0; i < log_fg_probs.size(); ++i) {
    log_fg_probs[i] = log(log_fg_probs[i]);
    log_bg_probs[i] = log(log_bg_probs[i]);
  }

  mixing = SplitDistro::log_sum_log_vec(log_fg_probs, log_fg_probs.size());
  const double bg_mixing =
    SplitDistro::log_sum_log_vec(log_bg_probs, log_bg_probs.size());
  const double mix_sum =
    ((mixing > bg_mixing) ? mixing + log(1 + exp(bg_mixing - mixing))
                          : bg_mixing + log(1 + exp(mixing - bg_mixing)));
  mixing = exp(mixing - mix_sum);
}

void
TwoStateSplitResolveMixture(
  const std::vector<double> &values, const std::vector<double> &vals_a,
  const std::vector<double> &vals_b, const std::vector<double> &scales,
  const std::size_t max_iterations, const double tolerance, const bool verbose,
  SplitDistro &fg_distro, SplitDistro &bg_distro, double &mixing) {
  // partition the observations to get the initial guess
  initialize_values(values, vals_a, vals_b, scales, fg_distro, bg_distro);

  std::vector<double> fg_probs(values.size(), 0);
  std::vector<double> bg_probs(values.size(), 0);
  mixing = 0.5;

  if (verbose)
    std::cout << '\n'
              << std::setw(10) << "DELTA" << std::setw(14) << "(PARAMS,MIX)\n";

  // Do the expectation maximization
  double prev_score = std::numeric_limits<double>::max();
  for (std::size_t itr = 0; itr < max_iterations; ++itr) {
    const double score = expectation_step(values, scales, mixing, fg_distro,
                                          bg_distro, fg_probs, bg_probs);
    maximization_step(/*values, */ vals_a, vals_b, scales, fg_probs, bg_probs,
                      mixing, fg_distro, bg_distro);

    if (verbose) {
      // clang-format off
      std::cout << std::setw(10) << std::setprecision(4)
                << (prev_score - score) / prev_score << "\t"
                << std::setw(14) << fg_distro.tostring() << " "
                << std::setw(10) << mixing << " "
                << std::setw(14) << bg_distro.tostring() << " "
                << std::setw(10) << 1 - mixing << " "
                << '\n';
      // clang-format on
    }
    if ((prev_score - score) / prev_score < tolerance)
      break;
    prev_score = score;
  }
}
