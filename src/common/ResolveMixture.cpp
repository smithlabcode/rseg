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

#include "ResolveMixture.hpp"

#include <memory>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <limits>

using std::vector;
using std::auto_ptr;
using std::cerr;
using std::endl;
using std::accumulate;


void
partition_values(const vector<double> &values, 
		 const size_t n_distros, 
		 vector<vector<double> > &partitioned) {

  const double val_weight = accumulate(values.begin(), values.end(), 0.0);
  assert(val_weight > 0);
  double val_per_part = std::ceil(val_weight/n_distros);
  
  vector<double> helper(values);
  sort(helper.begin(), helper.end(), std::greater<double>());
  
  partitioned.push_back(vector<double>());
  double sum = 0;
  for (size_t i = 0; i < helper.size(); ++i) {
    if (sum >= val_per_part) {
      partitioned.push_back(vector<double>());
      sum = 0;
    }
    partitioned.back().push_back(helper[i]);
    sum += helper[i];
  }
}


double
expectation_step(const vector<double> &values, 
		 const vector<double> &mixing,
		 const vector<Distro> &distros,
		 vector<vector<double> > &probs) {
  
  vector<double>::const_iterator x_idx = values.begin();
  const vector<double>::const_iterator x_lim = values.end();
  vector<vector<double> >::iterator y_idx = probs.begin();
  
  double score = 0;
  
  const size_t n_distros = distros.size();

  vector<double> parts(n_distros);
  vector<double> log_mixing(n_distros);
  for (size_t i = 0; i < n_distros; ++i) {
    log_mixing[i] = log(mixing[i]);
    assert(finite(log_mixing[i]));
  }
  for (size_t i = 0; i < values.size(); ++i) {
    
    copy(log_mixing.begin(), log_mixing.end(), parts.begin());
    for (size_t j = 0; j < n_distros; ++j) {
      parts[j] += distros[j].log_likelihood(values[i]);
      assert(finite(parts[j]));
    }
    
    const double denom = Distro::log_sum_log_vec(parts, parts.size());
    assert(finite(denom));
    for (size_t j = 0; j < n_distros; ++j)
      probs[j][i] = exp(parts[j] - denom);
    score += denom;
  }
  
  return score;
}


void
maximization_step(const vector<double> &values, 
		  const vector<vector<double> > &probs,
		  vector<double> &mixing,
		  vector<Distro> &distros) {
  const size_t n_distros = distros.size();
  const size_t vals_size = values.size();

  for (size_t i = 0; i < n_distros; ++i) {
    distros[i].estimate_params_ml(values, probs[i]);
    mixing[i] = Distro::log_sum_log_vec(probs[i], vals_size);
  }
  const double mix_sum = Distro::log_sum_log_vec(mixing, n_distros);
  for (vector<double>::iterator i(mixing.begin()); i != mixing.end(); ++i)
    *i = exp(*i - mix_sum);
}


void
ResolveMixture(const vector<double> &values,
	       const size_t max_iterations, const double tolerance, 
	       const int VERBOSE,
	       vector<Distro> &distros, vector<double> &mixing) {
  
  const size_t n_distros = distros.size();
  
  // partition the observations to get the initial guess
  vector<vector<double> > partitioned;
  partition_values(values, n_distros, partitioned);

  // Obtain initial estimates of distribution values
  for (size_t i = 0; i < n_distros; ++i)
    distros[i].estimate_params_ml(partitioned[i]);
  
  vector<vector<double> > probs(n_distros,
				vector<double>(values.size(), 0));
  mixing = vector<double>(n_distros, 1.0/n_distros);

  if (VERBOSE)
    cerr << endl << std::setw(10) << "DELTA"
	 << std::setw(14) << "(PARAMS,MIX)" << endl;
  
  // Do the expectation maximization
  double prev_score = std::numeric_limits<double>::max();
  for (size_t itr = 0; itr < max_iterations; ++itr) {
    const double score = expectation_step(values, mixing,
					  distros, probs);
    
    maximization_step(values, probs, mixing, distros);
    
    if (VERBOSE) {
      cerr << std::setw(10) << std::setprecision(4) 
	   << (prev_score - score)/prev_score << "\n";
      for (size_t i = 0; i < n_distros; ++i)
	cerr << "\t" 
	     << std::setw(14) << distros[i].tostring() << " " 
	     << std::setw(10) << mixing[i] << endl;
      cerr << endl;
    }
    if ((prev_score - score)/prev_score < tolerance)
      break;
    prev_score = score;
  }
}
