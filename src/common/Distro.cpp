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

#include "Distro.hpp"
#include "smithlab_utils.hpp"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <iostream>
#include <limits>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_diff.h>
#include <gsl/gsl_statistics_double.h>

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::bind2nd;
using std::divides;
using std::numeric_limits;
using std::pair;
using std::make_pair;
using std::accumulate;
using std::max_element;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// DISTRO__

double
Distro_::log_sum_log_vec(const vector<double> &vals, size_t limit) {
  const vector<double>::const_iterator x = 
    max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(finite(sum));
#endif
    }
  }
  return max_val + log(sum);
}

// Distro_::Distro_(const Distro_ &other) : 
//   params(other.params),
//   workspace_vals(other.workspace_vals),
//   workspace_probs(other.workspace_probs) {
//   rng = gsl_rng_clone(other.rng);
// }

Distro_::Distro_() {
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);
}

Distro_::Distro_(const std::vector<double> p) : params(p) {
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);
}

Distro_::~Distro_() {gsl_rng_free(rng);}
void 
Distro_::seed(int s) {gsl_rng_set(rng, s);}

string
Distro_::tostring() const {
  std::ostringstream os;
  if (!params.empty())
    os << std::setprecision(4) << params.front();
  for (size_t i = 1; i < params.size(); ++i)
    os << " " << std::setprecision(4) << params[i];
  return os.str();
}


double
Distro_::log_likelihood(std::vector<double>::const_iterator a,
			std::vector<double>::const_iterator b) const {
  double l = 0;
  for (; a < b; ++a)
    l += this->log_likelihood(*a);
  return l;
}

double 
Distro_::operator()(const double val) const {
  return exp(log_likelihood(val));
}


double 
Distro_::operator()(const vector<double> &vals) const {
  const size_t lim = vals.size();
  double l = 1;
  for (size_t i = 0; i < lim; ++i)
    l *= this->operator()(vals[i]);
  return l;
}


double 
Distro_::log_likelihood(const std::vector<double> &vals) const {
  double l = 0;
  const size_t lim = vals.size();
  for (size_t i = 0; i < lim; ++i)
    l += this->log_likelihood(vals[i]);
  return l;
}

double 
Distro_::log_likelihood(const std::vector<double> &vals,
                        const std::vector<double> &scales) const {
  double l = 0;
  const size_t lim = vals.size();
  for (size_t i = 0; i < lim; ++i)
      l += this->log_likelihood(vals[i], scales[i]);
  return l;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// DISTRO

bool
Distro::has_params(const string &name) {
  return (smithlab::split(name, ",").size() >= 2);
}

Distro::Distro(const std::string &n, const std::string &params) :
  name(n), d(distro_factory(n, params)) {}

Distro::Distro(const std::string &n, const std::vector<double> &params) :
  name(n), d(distro_factory(name)) {d->set_params(params);}

Distro::Distro(const std::string &s) : d(distro_factory(s)) {
    vector<string> name_split;
    if (s.find(",") == string::npos)  // whitespaces seperated 
        name_split = smithlab::split_whitespace_quoted(s);
    else // comma seperated
        name_split = smithlab::split(s, ",");
    name = name_split[0];
  name = name_split[0];
}

Distro::Distro(const Distro &rhs) : name(rhs.name), d(distro_factory(rhs.name)) {
  vector<double> tmp_params(rhs.get_params());
  d->set_params(tmp_params);
}

Distro& 
Distro::operator=(const Distro &rhs) {
  if (this != &rhs) {
    name = rhs.name;
    d = distro_factory(rhs.name);
    vector<double> tmp_params(rhs.get_params());
    d->set_params(tmp_params);
  }
  return *this;
}

Distro::~Distro() {if (d) delete d;}

double
Distro::operator()(double val) const {return (*d)(val);}

double
Distro::operator()(const std::vector<double> &vals) const {
  return (*d)(vals);
}

void
Distro::estimate_params_ml(const std::vector<double> &vals) {
  d->estimate_params_ml(vals);
}

void
Distro::estimate_params_ml(const std::vector<double> &vals,
                           const std::vector<double> &scales,
                           const std::vector<double> &probs) {
    d->estimate_params_ml(vals, scales, probs);
}

void
Distro::estimate_params_ml(const std::vector<double> &vals,
			   const std::vector<double> &weights) {
  d->estimate_params_ml(vals, weights);
}  

double
Distro::log_likelihood(const std::vector<double> &vals) const {
  return d->log_likelihood(vals);
}

double
Distro::log_likelihood(const std::vector<double> &vals,
                       const std::vector<double> &scales) const {
    return d->log_likelihood(vals, scales);
}



double
Distro::log_likelihood(std::vector<double>::const_iterator a,
		       std::vector<double>::const_iterator b) const {
  return d->log_likelihood(a, b);
}

double
Distro::log_likelihood(double val) const {return d->log_likelihood(val);}

double
Distro::log_likelihood(const double &val,
                       const double &scale) const
{return d->log_likelihood(val, scale);}


std::string
Distro::tostring() const {return name + string(" ") + d->tostring();}

std::ostream&
operator<<(std::ostream& s, const Distro& distro) {
  return s << distro.tostring();
}

double
Distro::log_sum_log_vec(const vector<double> &vals, size_t limit) {
  return Distro_::log_sum_log_vec(vals, limit);
}

////////////////////////////////////////////////////////////////////////
//
// DISTRO FACTORY


Distro_ *
distro_factory(string name, string params) {
  Distro_ *distro;
  if (name == "exp")
    distro = new ExpDistro();
  else if (name == "pois")
    distro = new PoisDistro();
  else if (name == "nbd")
    distro = new NegBinomDistro();
  else if (name == "geo")
    distro = new GeoDistro();
  else if (name == "emp")
    distro = new EmpiricalDistro();
  else if (name == "discemp")
    distro = new DiscEmpDistro();
  else throw SMITHLABException("bad distribution name \"" + name + "\"");
  
  vector<string> params_split = smithlab::split(params, ",");
  if (params_split.size() != distro->required_params())
    throw SMITHLABException("bad number of params: " + smithlab::toa(params_split.size()) +
			" for distro " + name + "\"");
  else {
    vector<double> params_vec;
    for (size_t i = 0; i < params_split.size(); ++i)
      params_vec.push_back(atof(params_split[i].c_str()));
    distro->set_params(params_vec);
  }
  return distro;
}


Distro_ *
distro_factory(string name_arg) {
    vector<string> name_split;
    if (name_arg.find(",") == string::npos)  // whitespaces seperated 
        name_split = smithlab::split_whitespace_quoted(name_arg);
    else // comma seperated
        name_split = smithlab::split(name_arg, ",");

  const string name = name_split.front();
  
  Distro_ *distro;
  if (name == "exp")
    distro = new ExpDistro();
  else if (name == "pois")
    distro = new PoisDistro();
  else if (name == "nbd")
    distro = new NegBinomDistro();
  else if (name == "geo")
    distro = new GeoDistro();
  else if (name == "emp")
    distro = new EmpiricalDistro();
  else if (name == "discemp")
    distro = new DiscEmpDistro();
  else throw SMITHLABException("bad distribution name \"" + name + "\"");

  if (name_split.size() > 1) {
    
    vector<string> params_split(vector<string>(name_split.begin() + 1,
					       name_split.end()));
    if (params_split.size() != distro->required_params())
      throw SMITHLABException("bad number of params: " + smithlab::toa(params_split.size()) +
			  " for distro " + name + "\"");
    else {
      vector<double> params_vec;
      for (size_t i = 0; i < params_split.size(); ++i)
	params_vec.push_back(atof(params_split[i].c_str()));
      distro->set_params(params_vec);
    }
  }
  return distro;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


double
ExpDistro::sample() const {
  assert(!params.empty());
  return gsl_ran_exponential(Distro_::rng, params.front());
}


double 
ExpDistro::log_likelihood(const double val) const {
  return -log(params[0]) - val/params[0];
}

double
ExpDistro::log_likelihood(const double &val,
                          const double &scale) const
{
    //// TEST NEEDED
    return -log(params[0] * scale) - val / params[0] / scale;
}


ExpDistro::ExpDistro(const ExpDistro &rhs) : 
  Distro_(rhs.params) {}

ExpDistro&
ExpDistro::operator=(const ExpDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
    rng = gsl_rng_alloc(gsl_rng_default);
  }
  return *this;
}


void
ExpDistro::estimate_params_ml(const vector<double> &vals) {
  params.front() = accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
}

void
ExpDistro::estimate_params_ml(const vector<double> &vals,
                              const vector<double> &scales,
                              const vector<double> &probs) {
    cerr << "ExpDistro::estimate_params_ml:  "
         << "Using an untested function" << endl;

    vector<double> values(vals.begin(), vals.end());
    for (size_t i = 0; i < values.size(); ++i)
        values[i] /= scales[i];
    if (probs.size() == 0)
        estimate_params_ml(vals);
    else
        estimate_params_ml(vals, probs);
}


void
ExpDistro::estimate_params_ml(const vector<double> &vals,
			      const vector<double> &probs) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }
  const double prob_sum = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(log_sum_log_vec(workspace_vals, lim))/prob_sum;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


double
PoisDistro::sample() const {
  assert(!params.empty());
  return gsl_ran_poisson(Distro_::rng, params.front());
}


double 
PoisDistro::log_likelihood(const double val) const {
  return -params.front() + val*log(params.front()) - 
    gsl_sf_lnfact(static_cast<size_t>(val));
}

double
PoisDistro::log_likelihood(const double &val,
                          const double &scale) const
{
    const double lambda = params[0] * scale;
    return -lambda + val*log(lambda) - 
        gsl_sf_lnfact(static_cast<size_t>(val));
}




PoisDistro::PoisDistro(const PoisDistro &rhs) : Distro_(rhs.params) {}


PoisDistro&
PoisDistro::operator=(const PoisDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
    rng = gsl_rng_alloc(gsl_rng_default);
  }
  return *this;
}


void
PoisDistro::estimate_params_ml(const vector<double> &vals) {
  params.front() = accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
}


void
PoisDistro::estimate_params_ml(const vector<double> &vals,
			       const vector<double> &probs) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }
  const double prob_sum = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(log_sum_log_vec(workspace_vals, lim))/prob_sum;
}

void
PoisDistro::estimate_params_ml(const vector<double> &vals,
                               const vector<double> &scales,
                               const vector<double> &probs) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i)
  {
      workspace_probs[i] = log(probs[i]) + log(scales[i]);
      workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }
  const double prob_sum = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(log_sum_log_vec(workspace_vals, lim))/prob_sum;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

const double
NegBinomDistro::max_allowed_alpha = 100;

const double
NegBinomDistro::min_allowed_alpha = 1e-20;

const double 
NegBinomDistro::alpha_allowed_error = 1e-10; 


void
NegBinomDistro::set_helpers() {
  n_helper = 1/params[1];
  p_helper = n_helper/(n_helper + params[0]);
  n_log_p_minus_lngamma_n_helper = n_helper*log(p_helper) - 
    gsl_sf_lngamma(n_helper);
  log_q_helper = log(1 - p_helper);
  //TODO: should check that these are valid!!!
}


void
NegBinomDistro::set_params(const std::vector<double> &p) {
  Distro_::set_params(p); 
  set_helpers();
}


double
NegBinomDistro::sample() const {
  assert(params.size() == 2);
  const double mu = params.front();
  const double one_over_alpha = 1/params.back();
  return gsl_ran_negative_binomial(Distro_::rng, 
				   one_over_alpha/(one_over_alpha + mu), 
				   one_over_alpha);
}


double 
NegBinomDistro::log_likelihood(const double val) const {
  const double P = (gsl_sf_lngamma(val + n_helper) - 
		    gsl_sf_lnfact(static_cast<size_t>(val))) +
    n_log_p_minus_lngamma_n_helper + val*log_q_helper;
  if (!finite(P))
    return -40;
  return P;
}

double
NegBinomDistro::log_likelihood(const double &val,
                               const double &scale) const
{
    //// TEST NEEDED
    // alpha is not scaled
    const double scaled_p_helper = n_helper / (n_helper + params[0] * scale);
    const double scaled_n_log_p_minus_lngamma_n_helper =
        n_helper*log(scaled_p_helper) - gsl_sf_lngamma(n_helper);
    const double scaled_log_q_helper = log(1 - scaled_p_helper);
    
    
    const double P = (gsl_sf_lngamma(val + n_helper) - 
                      gsl_sf_lnfact(static_cast<size_t>(val))) +
        scaled_n_log_p_minus_lngamma_n_helper + val * scaled_log_q_helper;
    if (!finite(P))
        return -40;
    return P;
}



NegBinomDistro::NegBinomDistro(const NegBinomDistro &rhs) : 
  Distro_(rhs.params) {
  set_helpers();
}


NegBinomDistro&
NegBinomDistro::operator=(const NegBinomDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
    rng = gsl_rng_alloc(gsl_rng_default);
    set_helpers();
  }
  return *this;
}


static inline double
score_fun_first_term(const vector<double> &vals_hist, 
		     const double mu, const double alpha) {
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
	inner_sum += j/(1 + alpha*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum;
}


static inline double
alpha_score_function(const vector<double> &vals_hist, const double mu, 
		     const double alpha, const double vals_count) {
  const double one_plus_alpha_mu = 1 + alpha*mu;
  return (score_fun_first_term(vals_hist, mu, alpha)/vals_count + 
	  (log(one_plus_alpha_mu)/alpha - mu)/alpha);
}


void
NegBinomDistro::estimate_params_ml(const vector<double> &vals) {
  // This is the mu
  params.front() = std::accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
  
  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.end());
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < vals.size(); ++i)
    ++vals_hist[static_cast<size_t>(vals[i])];
  
  const double vals_count = vals.size();
  
  const double mu = params.front();
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  
  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error && fabs((a_high - a_low)/max(a_high, a_low)) > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0) 
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  params[1] = a_mid;

  set_helpers();
//   const size_t lim = vals.size();

//   // This is the mu
//   params.front() = std::accumulate(vals.begin(), vals.begin() + lim, 0.0)/lim;
//   const double mu = params.front();
//   const double var = gsl_stats_variance_m(&vals.front(), 1, vals.size(), mu);
//   const double r = (mu*mu)/(var - mu);
//   // const double p = r/(r + params[0]);
//   params[1] = max(0.01, 1/r);
  
//   set_helpers();
}

void
NegBinomDistro::estimate_params_ml(const vector<double> &vals,
				   const vector<double> &probs) {
//   const size_t lim = vals.size();
//   if (workspace_vals.size() < lim) {
//     workspace_vals.resize(lim);
//     workspace_probs.resize(lim);
//   }
//   for (size_t i = 0; i < lim; ++i) {
//     workspace_probs[i] = log(probs[i]);
//     workspace_vals[i] = log(vals[i]) + log(probs[i]);
//   }
//   const double vals_count = exp(log_sum_log_vec(workspace_probs, lim));
//   const double mu = exp(log_sum_log_vec(workspace_vals, lim))/vals_count;
//   const double var = gsl_stats_wvariance_m(&vals.front(), 1, &probs.front(), 1, 
// 					   vals.size(), mu);
//   const double r = (mu*mu)/(var - mu);
//   // const double p = r/(r + params[0]);
//   params[0] = mu;
//   params[1] = max(0.01, 1/r);
//   set_helpers();
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    // assert(finite(probs[i]));
    workspace_probs[i] = log(probs[i]);// - centering_value;
    // assert(finite(workspace_probs[i]));
    workspace_vals[i] = log(vals[i]) + log(probs[i]);// - centering_value;
  }
  
  const double vals_count = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(log_sum_log_vec(workspace_vals, lim))/vals_count;
  
  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.begin() + lim);
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < lim; ++i)
    vals_hist[static_cast<size_t>(vals[i])] += probs[i];
  
  const double mu = params.front();
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  
  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error && 
	 fabs((a_high - a_low)/max(a_high, a_low)) > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    //     cerr << diff << "\t"
    // 	 << prev_val << "\t"
    // 	 << mid_val << "\t"
    // 	 << a_low << "\t"
    // 	 << a_mid << "\t"
    // 	 << a_high << endl;
    diff = std::fabs((prev_val - mid_val)/std::max(mid_val, prev_val));
    prev_val = mid_val;
  }
  params[1] = a_mid;

  set_helpers();
}

double
llh_deriative_rt_alpha(const vector<double> &vals,
                       const vector<double> &scales,
                       const vector<double> &probs,
                       const vector<double> &vals_hist,
                       const double mu,
                       const double alpha)
{
    const double first_term = score_fun_first_term(vals_hist, mu, alpha);
    
    const double mu_times_alpha = mu * alpha;
    const double alpha_inverse = 1 / alpha;
    const double alpha_square_inverse = pow(alpha_inverse, 2.0);
    
    double second_term = 0;
    for (size_t i = 0; i < vals.size(); ++i)
    {
        const double one_plus_extra = 1 + scales[i] * mu_times_alpha;
        second_term += probs[i] *
            (alpha_square_inverse * log(one_plus_extra) -
             scales[i] * mu * (alpha_inverse + vals[i]) / one_plus_extra);
    }

    return first_term + second_term;
}

void 
NegBinomDistro::estimate_params_ml(
    const std::vector<double> &vals,
    const std::vector<double> &scales,
    const std::vector<double> &probs)
{
    const size_t lim = vals.size();
    if (workspace_vals.size() < lim)
    {
        workspace_vals.resize(lim);
        workspace_probs.resize(lim);
    }
    for (size_t i = 0; i < lim; ++i)
    {
        workspace_probs[i] = log(probs[i]) + log(scales[i]);
        workspace_vals[i] = log(vals[i]) + log(probs[i]);
    }

    // this is mu
    const double mu = exp(log_sum_log_vec(workspace_vals, lim)) /
        exp(log_sum_log_vec(workspace_probs, lim));

    // Now for the alpha
    const double max_value = *std::max_element(vals.begin(), vals.end());
    vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
    for (size_t i = 0; i < vals.size(); ++i)
        vals_hist[static_cast<size_t>(vals[i])] += probs[i];
    
    double a_low = min_allowed_alpha;
    double a_high = max_allowed_alpha;
    
    double a_mid = max_allowed_alpha;
    double diff = std::numeric_limits<double>::max();
    double prev_val = std::numeric_limits<double>::max();
    while (diff > alpha_allowed_error &&
           fabs((a_high - a_low) / max(a_high, a_low)) > alpha_allowed_error)
    {
        a_mid = (a_low + a_high)/2;
        const double mid_val =
            llh_deriative_rt_alpha(vals, scales, probs, vals_hist, mu, a_mid);
        
        if (mid_val < 0) 
            a_high = a_mid;
        else
            a_low = a_mid;

        diff = std::fabs((prev_val - mid_val)/prev_val);
        prev_val = mid_val;
    }
    
    params[0] = mu;
    params[1] = a_mid;
    
    set_helpers();
}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


double
GeoDistro::sample() const {
  assert(!params.empty());
  return gsl_ran_geometric(Distro_::rng, params.front());
}


double 
GeoDistro::log_likelihood(const double val) const {
  return log(params[0]) + (val - 1)*log(1 - params[0]);
}

double
GeoDistro::log_likelihood(const double &val,
                          const double &scale) const
{
    cerr << "GeoDistro::log_likelihood:  "
         << "Using an untested function" << endl;

    return 0;
}



GeoDistro::GeoDistro(const GeoDistro &rhs) : 
  Distro_(rhs.params) {}

GeoDistro&
GeoDistro::operator=(const GeoDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
    rng = gsl_rng_alloc(gsl_rng_default);
  }
  return *this;
}


void
GeoDistro::estimate_params_ml(const vector<double> &vals) {
  params.front() = vals.size()/accumulate(vals.begin(), vals.end(), 0.0);
  assert(params.front() <= 1.0);
}


void
GeoDistro::estimate_params_ml(const vector<double> &vals,
			      const vector<double> &probs) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }
  const double prob_sum = exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = prob_sum/exp(log_sum_log_vec(workspace_vals, lim));
}

void 
GeoDistro::estimate_params_ml(
    const std::vector<double> &vals,
    const std::vector<double> &scales,
    const std::vector<double> &probs)
{
    cerr << "GeoDistro::estimate_params_ml:  "
         << "Using an untested function" << endl;
    vector<double> values(vals);
    for (size_t i = 0; i < values.size(); ++i)
        values[i] /= scales[i];
    estimate_params_ml(values, probs);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include "Smoothing.hpp"

const double
EmpiricalDistro::MIN_PROB = 1e-10;


double
EmpiricalDistro::estimate_number_of_classes(const vector<double> &vals) {
  return min(sqrt(vals.size()), *max_element(vals.begin(), vals.end()) + 1);
}


double
EmpiricalDistro::estimate_bandwidth(const vector<double> &vals) {
  const size_t lim = vals.size();
  return gsl_stats_sd(&vals.front(), 1, lim);
}


double
EmpiricalDistro::estimate_number_of_classes(const vector<double> &vals,
					    const vector<double> &probs) {
  return min(sqrt(vals.size()), *max_element(vals.begin(), vals.end()) + 1);
}


double
EmpiricalDistro::estimate_bandwidth(const vector<double> &vals,
				    const vector<double> &probs) {
  
  assert(vals.size() == probs.size());
  
  const size_t lim = vals.size();
  
  // bandwidth parameter is the (weighted) standard deviation
  double x_tilda = 0;
  double xx_tilda = 0;
  for (size_t i = 0; i < lim; ++i) {
    const double x = vals[i];
    const double w = probs[i];
    x_tilda += w*x;
    xx_tilda += w*x*x;
  }
  x_tilda /= lim;
  xx_tilda /= lim;
  return std::sqrt(xx_tilda - x_tilda*x_tilda);
}


size_t
EmpiricalDistro::find_bin(const vector<double> &bins,
			  const double val) {
  const size_t bin =
    upper_bound(bins.begin(), bins.end(), val) - bins.begin() - 1;
  assert(bin < bins.size() && bin >= 0);
  return bin;
}


void
EmpiricalDistro::make_cumulative(std::vector<double> &vals) {
  for (size_t i = 1; i < vals.size(); ++i)
    vals[i] += vals[i - 1];
}


void
EmpiricalDistro::get_breaks(vector<double> data,
			    size_t n_vals, size_t n_class, double max_val,
			    vector<double> &breaks) {
  
  sort(data.begin(), data.end());
  // the max number of observations per class
  const double max_obs = data.size()/n_class;
  // the minimum span of a non-full class
  const double min_span = max_val/n_class;
  
  breaks.reserve(n_class + 2);
  breaks.push_back(0);
  
  size_t data_id = 0;
  for (size_t i = 0; data_id < data.size(); ++i) {
    
    double obs_count = 0;
    bool changed = false;
    while (data_id < data.size() && 
	   (!changed || obs_count < max_obs) && 
	   (data[data_id] - breaks.back()) < min_span) {
      ++obs_count;
      ++data_id;
      if (data[data_id] != data[data_id - 1]) 
	changed = true;
    }
    if (obs_count >= max_obs)
      breaks.push_back(data[data_id]);
    else breaks.push_back(breaks.back() + min_span);
  }
}


void
EmpiricalDistro::get_breaks(const vector<double> &in_data,
			    const vector<double> &weights,
			    size_t n_vals,
			    size_t n_class,
			    double max_val,
			    vector<double> &breaks) {
  
  vector<pair<double, double> > data;
  for (size_t i = 0; i < in_data.size(); ++i)
    data.push_back(make_pair(in_data[i], weights[i]));
  
  sort(data.begin(), data.end());
  // the max number of observations per class
  const double max_obs = accumulate(weights.begin(), weights.end(), 0.0)/n_class;
  // the minimum span of a non-full class
  const double min_span = (max_val + 1)/n_class;
  breaks.reserve(n_class + 2);
  breaks.push_back(0);
  
  size_t data_id = 0;
  for (size_t i = 0; data_id < data.size(); ++i) {
    
    double obs_count = 0;
    bool changed = false;
    while (data_id < data.size() && 
	   (!changed || obs_count < max_obs) && 
	   (data[data_id].first - breaks.back()) < min_span) {
      obs_count += data[data_id].second;
      ++data_id;
      if (data[data_id].first != data[data_id - 1].first) 
	changed = true;
    }
    if (obs_count >= max_obs)
      breaks.push_back(data[data_id].first);
    else
      breaks.push_back(breaks.back() + min_span);
  }
}


void
EmpiricalDistro::make_hist(const vector<double> &data,
			   size_t n_vals,
			   size_t n_class,
			   double max_val,
			   vector<double> &breaks,
			   vector<double> &hist) {

  get_breaks(data, n_vals, n_class, max_val, breaks);
  // add the final break
  breaks.push_back(numeric_limits<double>::max());
  
  // fill the counts
  hist.clear();
  hist.resize(breaks.size() - 1);
  for (size_t i = 0; i < n_vals; ++i) {
    const size_t id = find_bin(breaks, data[i]);
    assert(id < hist.size() && data[i] >= breaks[id] && data[i] < breaks[id + 1]);
    hist[id]++;
  }
  hist.back() = 0;
}


void
EmpiricalDistro::make_weighted_hist(const vector<double> &data,
				    const vector<double> &weights,
				    size_t n_vals,
				    size_t n_class,
				    double max_val,
				    vector<double> &breaks,
				    vector<double> &hist) {

  breaks.clear();
  get_breaks(data, weights, n_vals, n_class, max_val, breaks);
  // add the final break
  breaks.push_back(numeric_limits<double>::max());
  
  // fill the counts
  hist.clear();
  hist.resize(breaks.size() - 1);
  for (size_t i = 0; i < n_vals; ++i) {
    const size_t id = find_bin(breaks, data[i]);
    assert(id < hist.size());
    assert(data[i] >= breaks[id]);
    assert(data[i] < breaks[id + 1]);
    hist[id] += weights[i];
  }
  hist.back() = 0;
}



double
EmpiricalDistro::sample() const {
  assert(!params.empty());
  const size_t id = find_bin(cumulative, gsl_rng_uniform(rng));
  return gsl_ran_flat(rng, breaks[id], breaks[id + 1]);
}


double 
EmpiricalDistro::log_likelihood(const double val) const {
  return log_hist[find_bin(breaks, val)];
}

double
EmpiricalDistro::log_likelihood(const double &val,
                              const double &scale) const
{
    cerr << "EmpiricalDistro::log_likelihood:  "
         << "Using an untested function" << endl;
    return log_hist[static_cast<size_t>(val / scale)];
}




EmpiricalDistro::EmpiricalDistro(const EmpiricalDistro &rhs) : 
  Distro_(rhs.params), 
  breaks(rhs.breaks), 
  log_hist(rhs.log_hist),
  hist(rhs.hist),
  cumulative(rhs.cumulative) {}


EmpiricalDistro&
EmpiricalDistro::operator=(const EmpiricalDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
    // copy tables
    breaks = rhs.breaks;
    hist = rhs.hist;
    log_hist = rhs.log_hist;
    cumulative = rhs.cumulative;
  }
  return *this;
}


void
EmpiricalDistro::estimate_params_ml(const vector<double> &vals) {

  const size_t lim = vals.size();
    
  params[0] = estimate_bandwidth(vals);
  params[1] = estimate_number_of_classes(vals);
  
  // build histogram
  make_hist(vals, lim, static_cast<size_t>(params[1]),
	    *max_element(vals.begin(), vals.begin() + lim),
	    breaks, hist);
  
  // make the mids
  vector<double> mids(hist.size());
  for (size_t i = 0; i < mids.size(); ++i)
    mids[i] = (breaks[i] + breaks[i + 1])/2;
  
  // smooth histogram
  vector<double> smooth_hist;
  LocalLinearRegression(params[0], mids, hist, mids, smooth_hist);
  hist.swap(smooth_hist);
  smooth_hist.clear();
  
  // normalize the table for probs
  const double total = accumulate(hist.begin(), hist.end(), 0.0);
  transform(hist.begin(), hist.end(), hist.begin(), 
	    bind2nd(divides<double>(), total));

  for (size_t i = 0; i < hist.size(); ++i)
    hist[i] = max(MIN_PROB, hist[i]);
  
  // preprocess cumulative table lookup
  cumulative = hist;
  make_cumulative(cumulative);
  
  // preprocess log prob table lookup
  log_hist.resize(hist.size());
  for (size_t i = 0; i < hist.size(); ++i)
    log_hist[i] = max(log(hist[i]) - log(breaks[i + 1] - breaks[i]), log(MIN_PROB));
}



void
EmpiricalDistro::estimate_params_ml(const vector<double> &vals,
				    const vector<double> &probs) {

  assert(vals.size() == probs.size());

  const size_t lim = vals.size();
  
  params[0] = estimate_bandwidth(vals, probs);
  params[1] = estimate_number_of_classes(vals, probs);
  
  // build histogram
  make_weighted_hist(vals, probs, lim, static_cast<size_t>(params[1]),
		     *max_element(vals.begin(), vals.begin() + lim),
		     breaks, hist);
  
  // make the mids
  vector<double> mids(hist.size());
  for (size_t i = 0; i < mids.size(); ++i)
    mids[i] = (breaks[i] + breaks[i + 1])/2;
  
  // smooth histogram
  vector<double> smooth_hist;
  LocalLinearRegression(params[0], mids, hist, mids, smooth_hist);
  hist.swap(smooth_hist);
  smooth_hist.clear();
  
  // normalize the table for probs
  const double total = accumulate(hist.begin(), hist.end(), 0.0);
  transform(hist.begin(), hist.end(), hist.begin(), 
	    bind2nd(divides<double>(), total));
  
  for (size_t i = 0; i < hist.size(); ++i)
    hist[i] = max(MIN_PROB, hist[i]);
  
  // preprocess cumulative table lookup
  cumulative = hist;
  make_cumulative(cumulative);
  
  // preprocess log prob table lookup
  log_hist.resize(hist.size());
  for (size_t i = 0; i < hist.size(); ++i)
    log_hist[i] = max(log(hist[i]) - log(breaks[i + 1] - breaks[i]), log(MIN_PROB));
}

void 
EmpiricalDistro::estimate_params_ml(
    const std::vector<double> &vals,
    const std::vector<double> &scales,
    const std::vector<double> &probs)
{
    cerr << "EmpiricalDistro::estimate_params_ml:  "
         << "Using an untested function" << endl;
    vector<double> values(vals);
    for (size_t i = 0; i < values.size(); ++i)
        values[i] /= scales[i];
    estimate_params_ml(values, probs);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////
////// DISCRETE EMPIRICAL DISTRIBUTION
//////

const double DiscEmpDistro::MIN_PROB = 1e-20;

size_t
DiscEmpDistro::find_bin(const vector<double> &bins, const double val) {
  const size_t bin = upper_bound(bins.begin(), bins.end(), val) - bins.begin() - 1;
  assert(bin < bins.size() && bin >= 0);
  return bin;
}

void
DiscEmpDistro::make_cumulative(std::vector<double> &vals) {
  for (size_t i = 1; i < vals.size(); ++i)
    vals[i] += vals[i - 1];
}

void
DiscEmpDistro::make_hist(const vector<double> &data, size_t n_vals, 
			 size_t n_classes, double max_val, vector<double> &hist) {
  // fill the counts
  hist.clear();
  hist.resize(n_classes);
  for (size_t i = 0; i < n_vals; ++i) {
    assert(data[i] >= 0 && data[i] <= max_val);
    hist[static_cast<size_t>(data[i])]++;
  }
  hist.back() = 0;
}

void
DiscEmpDistro::make_weighted_hist(const vector<double> &data,
				  const vector<double> &weights, size_t n_vals, 
				  size_t n_classes, double max_val, vector<double> &hist) {
  // fill the counts
  hist.clear();
  hist.resize(n_classes);
  for (size_t i = 0; i < n_vals; ++i) {
    assert(data[i] >= 0 && data[i] <= n_classes);
    hist[static_cast<size_t>(data[i])] += weights[i];
  }
  hist.back() = 0;
}

double
DiscEmpDistro::sample() const {
  assert(!params.empty());
  const size_t id = find_bin(cumulative, gsl_rng_uniform(rng));
  assert(id >= 0 && id <= max_val);
  return id;
}

double 
DiscEmpDistro::log_likelihood(const double val) const {
  return log_hist[static_cast<size_t>(val)];
}

double
DiscEmpDistro::log_likelihood(const double &val,
                              const double &scale) const
{
    cerr << "DiscEmpDistro::log_likelihood:  "
         << "Using an untested function" << endl;
    return log_hist[static_cast<size_t>(val / scale)];
}


DiscEmpDistro::DiscEmpDistro(const DiscEmpDistro &rhs) : 
  Distro_(rhs.params), 
  log_hist(rhs.log_hist),
  hist(rhs.hist),
  cumulative(rhs.cumulative) {}

DiscEmpDistro&
DiscEmpDistro::operator=(const DiscEmpDistro &rhs) {
  if (this != &rhs) {
    this->Distro_::params = rhs.Distro_::params;
    // copy tables
    hist = rhs.hist;
    log_hist = rhs.log_hist;
    cumulative = rhs.cumulative;
  }
  return *this;
}


void
DiscEmpDistro::estimate_params_ml(const vector<double> &vals) {

  const size_t lim = vals.size();

  params[0] = *max_element(vals.begin(), vals.begin() + lim);
  
  n_classes = static_cast<size_t>(params[0]) + 1;
  max_val = params[0];

  params[1] = accumulate(vals.begin(), vals.begin() + lim, 0.0)/lim;
  
  // build histogram
  make_hist(vals, lim, n_classes, max_val, hist);
  
  // make the mids
  vector<double> mids(hist.size());
  for (size_t i = 0; i < mids.size(); ++i)
    mids[i] = i;// + 0.5;
  
  // smooth histogram
  vector<double> smooth_hist;
  LocalLinearRegression(3.0, mids, hist, mids, smooth_hist);
  hist.swap(smooth_hist);
  smooth_hist.clear();
  
  // normalize the table for probs
  const double total = accumulate(hist.begin(), hist.end(), 0.0);
  transform(hist.begin(), hist.end(), hist.begin(), 
	    bind2nd(divides<double>(), total));
  
  for (size_t i = 0; i < hist.size(); ++i)
    hist[i] = max(MIN_PROB, hist[i]);
  
  // preprocess cumulative table lookup
  cumulative = hist;
  make_cumulative(cumulative);
  
  // preprocess log prob table lookup
  log_hist.resize(hist.size());
  for (size_t i = 0; i < hist.size(); ++i)
    log_hist[i] = log(hist[i]);
}

void
DiscEmpDistro::estimate_params_ml(const vector<double> &vals,
				  const vector<double> &probs) {

  assert(vals.size() == probs.size());

  const size_t lim = vals.size();
  
  params[0] = *max_element(vals.begin(), vals.begin() + lim);
  
  n_classes = static_cast<size_t>(params[0]) + 1;
  max_val = params[0];
  
  // build histogram
  make_weighted_hist(vals, probs, lim, n_classes, max_val, hist);

  params[1] = inner_product(vals.begin(), vals.begin() + lim,
			    probs.begin(), 0.0)/
    accumulate(probs.begin(), probs.begin() + lim, 0.0);
  
  // make the mids
  vector<double> mids(hist.size());
  for (size_t i = 0; i < mids.size(); ++i)
    mids[i] = i;// + 0.5;
  
  // smooth histogram
  vector<double> smooth_hist;
  LocalLinearRegression(1.5, mids, hist, mids, smooth_hist);
  hist.swap(smooth_hist);
  smooth_hist.clear();
  
  // normalize the table for probs
  const double total = accumulate(hist.begin(), hist.end(), 0.0);
  transform(hist.begin(), hist.end(), hist.begin(), 
	    bind2nd(divides<double>(), total));
  
  for (size_t i = 0; i < hist.size(); ++i)
    hist[i] = max(MIN_PROB, hist[i]);
  
  // preprocess cumulative table lookup
  cumulative = hist;
  make_cumulative(cumulative);
  
  // preprocess log prob table lookup
  log_hist.resize(hist.size());
  for (size_t i = 0; i < hist.size(); ++i)
    log_hist[i] = log(hist[i]);
}

void 
DiscEmpDistro::estimate_params_ml(
    const std::vector<double> &vals,
    const std::vector<double> &scales,
    const std::vector<double> &probs)
{
    cerr << "DiscEmpDistro::estimate_params_ml:  "
         << "Using an untested function" << endl;
    vector<double> values(vals);
    for (size_t i = 0; i < values.size(); ++i)
        values[i] /= scales[i];
    estimate_params_ml(values, probs);
}


  
//   // Now for the alpha
//   const double max_value = *max_element(vals.begin(), vals.begin() + lim);
//   vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
//   for (size_t i = 0; i < lim; ++i)
//     vals_hist[static_cast<size_t>(vals[i])] += probs[i];
  
//   double a_low = min_allowed_alpha;
//   double a_high = max_allowed_alpha;
  
//   double a_mid = max_allowed_alpha;
//   double diff = numeric_limits<double>::max();
//   double prev_val = numeric_limits<double>::max();
//   while (diff > alpha_allowed_error && fabs(a_high - a_low) > alpha_allowed_error) {
//     a_mid = (a_low + a_high)/2;
//     const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
//     if (mid_val < 0) a_high = a_mid;
//     else a_low = a_mid;
//     diff = std::fabs((prev_val - mid_val)/std::max(mid_val, prev_val));
//     prev_val = mid_val;
//   }
//   params[1] = a_mid;

  // Now for the alpha
//   const double max_value = *max_element(vals.begin(), vals.begin() + lim);
//   vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
//   for (size_t i = 0; i < lim; ++i)
//     ++vals_hist[static_cast<size_t>(vals[i])];
  
//   const double vals_count = lim;
  
//   const double mu = params.front();
//   double a_low = min_allowed_alpha;
//   double a_high = max_allowed_alpha;
  
//   double a_mid = max_allowed_alpha;
//   double diff = numeric_limits<double>::max();
//   double prev_val = numeric_limits<double>::max();
//   while (diff > alpha_allowed_error && fabs(a_high - a_low) > alpha_allowed_error) {
//     a_mid = (a_low + a_high)/2;
//     const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
//     if (mid_val < 0) a_high = a_mid;
//     else a_low = a_mid;
//     diff = std::fabs((prev_val - mid_val)/prev_val);
//     prev_val = mid_val;
//   }
//   params[1] = a_mid;


//   const double a = mu*mu;
//   const double b = 2*mu;
//   const double c = (1 - var);
//   const double discrim = b*b - 4*a*c;
//   const double solution = (-b + std::sqrt(discrim))/(2*a);
//   const double r = 1/solution;
//   params[1] = 1/r;
//   const double mu = params.front();
//   const double var = gsl_stats_variance(&vals.front(), 1, 
// 					vals.size());
  
//   const double r = 1.0/((-(2*mu) + std::sqrt((2*mu)*(2*mu) - 4*((mu*mu)*(1 - var))))/
//  			(2*mu*mu));
//   //   double discrim_root = fabs(std::sqrt((2*mu)*(2*mu) - 4*((mu*mu)*(1 - var))));
//   //   double r = 1.0/((-(2*mu) + discrim_root)/(2*mu*mu));

//   cerr << "p=" << p << "\t" << "r=" << r << "\n" 
//        << "mu=" << mu << "\t" << r*(1 - p)/p << "\n"
//        << "var=" << var << "\t" << r*(1 - p)/(p*p) << "\n";
  

  
//   set_helpers();
  
  //   double r = 1/params[1];
  //   double p = r/(r + params[0]);
  
  //   cerr << "p=" << p << "\n" 
  //        << "r=" << r << "\n" 
  //        << "mu=" << params[0] << "\t" << r*(1 - p)/p << "\n"
  //        << "var=" << variance << "\t" << r*(1 - p)/(p*p) << "\n"
  //        << "var=" << variance << "\t" << r*(1 - p)/(p*p) << "\n";
  
  //   double rr = 1.0/((-(2*mu) + std::sqrt((2*mu)*(2*mu) - 4*((mu*mu)*(1 - variance))))/
  // 		   (2*mu*mu));
  
  //   cerr << rr << "\t" << r*(1 - p)/(p*p) << endl;
  //   cerr << mu << "\t" << r*(1 - p)/p << endl;
  // exit(0);
