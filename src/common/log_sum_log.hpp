/* Copyright (C) 2025 Andrew D Smith
 * Author: Andrew D. Smith
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

#ifndef NUMERICAL_UTILS_HPP_
#define NUMERICAL_UTILS_HPP_

#include <algorithm>
#include <cmath>
#include <vector>

[[nodiscard]] inline double
log_sum_log(const double p, const double q) {
  if (p == 0)
    return q;
  if (q == 0)
    return p;
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + std::log(1.0 + std::exp(smaller - larger));
}

[[nodiscard]] inline double
log_sum_log(const double p, const double q, const double r) {
  return log_sum_log(log_sum_log(p, q), r);
}

[[nodiscard]] inline auto
log_sum_log_vec(const std::vector<double> &vals,
                const std::size_t limit) -> double {
  const auto x = std::max_element(std::cbegin(vals), std::cbegin(vals) + limit);
  const double max_val = *x;
  const std::size_t max_idx = x - std::cbegin(vals);
  double sum = 1.0;
  for (std::size_t i = 0; i < limit; ++i)
    if (i != max_idx)  // remove if vectorizing
      sum += std::exp(vals[i] - max_val);
  return max_val + std::log(sum);
}

[[nodiscard]] inline auto
log_sum_log(const std::vector<double>::const_iterator &begin,
            const std::vector<double>::const_iterator &end) -> double {
  const auto max_itr = std::max_element(begin, end);
  const double max_val = *max_itr;

  double sum = 1.0;
  for (auto itr = begin; itr != end; ++itr)
    if (itr != max_itr)
      sum += std::exp(*itr - max_val);

  return max_val + std::log(sum);
}

#endif  // NUMERICAL_UTILS_HPP_
