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

#ifndef NUMERICAL_UTILS_HPP
#define NUMERICAL_UTILS_HPP

#include "numerical_utils.cpp"

#include <cmath>
#include <vector>
#include <algorithm>

using std::vector;

double
log_sum_log(const double p, const double q) {
  if (p == 0) {
    return q;
  }
  else if (q == 0) {
    return p;
  }
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


double
log_sum_log(const double p, const double q, const double r) {
  if (p == 0) {
    return log_sum_log(q, r);
  }
  else if (q == 0) {
    return log_sum_log(p, r);
  }
  else if (r == 0) {
    return log_sum_log(p, r);
  }
  else {
    return log_sum_log(p, log_sum_log(q, r));
  }
}

double
log_sum_log_vec(const vector<double> &vals, const size_t limit) {
  const vector<double>::const_iterator x = 
    std::max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
    }
  }
  return max_val + log(sum);
}

#endif

