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

#include "SplitDistro.hpp"
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
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_multimin.h>

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
using std::ios_base;

static double
split_distro_log_sum_log_vec(const vector<double> &vals, size_t limit) {
  const vector<double>::const_iterator x = 
    std::max_element(vals.begin(), vals.begin() + limit);
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

double
SplitDistro_::log_likelihood(std::vector<double>::const_iterator a,
			     std::vector<double>::const_iterator b) const {
  double l = 0;
  for (; a < b; ++a)
    l += this->log_likelihood(*a);
  return l;
}

double 
SplitDistro_::operator()(const double val) const {
  return exp(log_likelihood(val));
}

double 
SplitDistro_::operator()(const vector<double> &vals) const {
  const size_t lim = vals.size();
  double l = 1;
  for (size_t i = 0; i < lim; ++i)
    l *= this->operator()(vals[i]);
  return l;
}

double 
SplitDistro_::log_likelihood(const std::vector<double> &vals) const {
  double l = 0;
  const size_t lim = vals.size();
  for (size_t i = 0; i < lim; ++i)
    l += this->log_likelihood(vals[i]);
  return l;
}

double 
SplitDistro_::log_likelihood(const std::vector<double> &vals,
                             const std::vector<double> &scales) const {
  double l = 0;
  const size_t lim = vals.size();
  for (size_t i = 0; i < lim; ++i)
      l += this->log_likelihood(vals[i], scales[i]);
  return l;
}

SplitDistro::SplitDistro(const std::string &n, const std::string &params) :
  name(n), d(split_distro_factory(n, params)) {}

SplitDistro::SplitDistro(const std::string &n, const std::vector<double> &params) :
  name(n), d(split_distro_factory(name)) {d->set_params(params);}

SplitDistro::SplitDistro(const std::string &s) : d(split_distro_factory(s)) {
    vector<string> name_split;
    if (s.find(",") == string::npos)  // whitespaces seperated 
        name_split = smithlab::split_whitespace_quoted(s);
    else // comma seperated
        name_split = smithlab::split(s, ",");
    name = name_split[0];
}

SplitDistro::SplitDistro(const SplitDistro &rhs) : name(rhs.name), 
						   d(split_distro_factory(rhs.name)) {
  vector<double> tmp_params(rhs.get_params());
  d->set_params(tmp_params);
}

SplitDistro& 
SplitDistro::operator=(const SplitDistro &rhs) {
  if (this != &rhs) {
    name = rhs.name;
    d = split_distro_factory(rhs.name);
    vector<double> tmp_params(rhs.get_params());
    d->set_params(tmp_params);
  }
  return *this;
}

SplitDistro::~SplitDistro() {if (d) delete d;}

double
SplitDistro::operator()(double val) const {return (*d)(val);}

double
SplitDistro::operator()(const std::vector<double> &vals) const {
  return (*d)(vals);
}


void
SplitDistro::estimate_params_ml(const std::vector<double> &vals_a, 
				const std::vector<double> &vals_b) {
  d->estimate_params_ml(vals_a, vals_b);
}


void
SplitDistro::estimate_params_ml(const std::vector<double> &vals_a,
				const std::vector<double> &vals_b,
				const std::vector<double> &weights) {
  d->estimate_params_ml(vals_a, vals_b, weights);
}  

void
SplitDistro::estimate_params_ml(const std::vector<double> &vals_a,
                                const std::vector<double> &vals_b,
                                const std::vector<double> &scales,
                                const std::vector<double> &weights) {
    d->estimate_params_ml(vals_a, vals_b, scales, weights);
}  


double
SplitDistro::log_likelihood(const std::vector<double> &vals) const {
  return d->log_likelihood(vals);
}

double
SplitDistro::log_likelihood(const std::vector<double> &vals,
                            const std::vector<double> &scales) const {
    return d->log_likelihood(vals, scales);
}


double
SplitDistro::log_likelihood(std::vector<double>::const_iterator a,
			    std::vector<double>::const_iterator b) const {
  return d->log_likelihood(a, b);
}


double
SplitDistro::log_likelihood(double val) const {
  return d->log_likelihood(val);
}

double
SplitDistro::log_likelihood(const double val, const double scale) const {
    return d->log_likelihood(val, scale);
}

std::string
SplitDistro::tostring() const {
  return name + string(" ") + d->tostring();
}



std::ostream&
operator<<(std::ostream& s, const SplitDistro& distro) {
  return s << distro.tostring();
}


double
SplitDistro::log_sum_log_vec(const vector<double> &vals, size_t limit) {
  return split_distro_log_sum_log_vec(vals, limit);
}

////////////////////////////////////////////////////////////////////////
//
// DISTRO FACTORY


SplitDistro_ *
split_distro_factory(string name, string params) {
  SplitDistro_ *distro;
  if (name == "skel")
    distro = new SkellamDistro();
  else if (name == "gauss")
    distro = new GaussianDistro();
  else if (name == "nbdiff")
    distro = new NegBinomDiffDistro();
  else if (name == "discdiff")
    distro = new DiscEmpSplitDistro();
  else {
    cerr << "bad distribution name \"" << name << "\"" << endl;
    exit(EXIT_FAILURE);
  }
  
  vector<string> params_split = smithlab::split(params, ",");
  if (params_split.size() != distro->required_params())
    throw SMITHLABException("bad number of params: " + smithlab::toa(params_split.size()) +
			" for split distro " + name + "\"");
  else {
    vector<double> params_vec;
    for (size_t i = 0; i < params_split.size(); ++i) {
      params_vec.push_back(atof(params_split[i].c_str()));
    }
    distro->set_params(params_vec);
  }
  return distro;
}


SplitDistro_ *
split_distro_factory(string name_arg) {

    vector<string> name_split;
    if (name_arg.find(",") == string::npos)  // whitespaces seperated 
        name_split = smithlab::split_whitespace_quoted(name_arg);
    else // comma seperated
        name_split = smithlab::split(name_arg, ",");
    
  const string name = name_split.front();
  
  SplitDistro_ *distro;
  if (name == "skel")
    distro = new SkellamDistro();
  else if (name == "nbdiff")
    distro = new NegBinomDiffDistro();
  else if (name == "gauss")
    distro = new GaussianDistro();
  else if (name == "discdiff")
    distro = new DiscEmpSplitDistro();
  else
    throw SMITHLABException("bad split distribution name \"" + name + "\"");

  if (name_split.size() > 1) {
    vector<string> params_split(vector<string>(name_split.begin() + 1,
					       name_split.end()));
    if (params_split.size() != distro->required_params())
      throw SMITHLABException("bad number of params: " + smithlab::toa(params_split.size()) +
			  " for split distro " + name + "\"");
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
//
// DISTRO

SplitDistro_::SplitDistro_() {
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);
}

SplitDistro_::SplitDistro_(const std::vector<double> p) : params(p) {
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);
}

SplitDistro_::~SplitDistro_() {
  gsl_rng_env_setup();
  gsl_rng_free(rng);
}

void
SplitDistro_::seed(int s) {
  gsl_rng_set(rng, s);
}

string
SplitDistro_::tostring() const {
  const ios_base::fmtflags split_distro_opts = ios_base::showpoint;
  std::ostringstream os;
  os.flags(split_distro_opts);
  if (!params.empty())
    os << std::setprecision(3) << params.front();
  for (size_t i = 1; i < params.size(); ++i)
    os << " " << std::setprecision(3) << params[i];
  return os.str();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

double
SkellamDistro::sample() const {
  assert(params.size() == 2);
  return gsl_ran_poisson(SplitDistro_::rng, params.front()) -
    gsl_ran_poisson(SplitDistro_::rng, params.back());
}

double 
SkellamDistro::log_likelihood(const double val) const {
  if (std::fabs(val) > 120)
    return -40;

//   const double r = -(params.front() + params.back()) + 
//     (val/2.0)*(log(params.front()) - log(params.back())) +
//     log(gsl_sf_bessel_In(static_cast<int>(val), 
// 			 2*sqrt(params.front()*params.back())));

  const double r = -(params.front() + params.back()) + 
      (val/2.0)*(log(params.front()) - log(params.back())) +
      log(gsl_sf_bessel_In(static_cast<int>(fabs(val)), 
                           2*sqrt(params.front()*params.back())));

  return r;
}

double 
SkellamDistro::log_likelihood(const double val,
                              const double scale) const {
  if (std::fabs(val) > 120)
    return -40;

  const double lambda1 = params[0] * scale;
  const double lambda2 = params[1] * scale;
  
  const double r = -(lambda1 + lambda2) + 
      (val/2.0)*(log(lambda1) - log(lambda2)) +
      log(gsl_sf_bessel_In(static_cast<int>(fabs(val)), 
                           2*sqrt(lambda1 * lambda2)));
  return r;
}

SkellamDistro::SkellamDistro(const SkellamDistro &rhs) : 
  SplitDistro_(rhs.params) {}

SkellamDistro&
SkellamDistro::operator=(const SkellamDistro &rhs) {
  if (this != &rhs) {
    this->SplitDistro_::params = rhs.SplitDistro_::params;
    rng = gsl_rng_alloc(gsl_rng_default);
  }
  return *this;
}

void
SkellamDistro::estimate_params_ml(const vector<double> &vals_a, 
				  const vector<double> &vals_b) {
  assert(vals_a.size() == vals_b.size());
  const size_t lim = vals_a.size();
  params.front() = accumulate(vals_a.begin(), vals_a.end(), 0.0)/lim;
  params.back() = accumulate(vals_b.begin(), vals_b.end(), 0.0)/lim;
}

void
SkellamDistro::estimate_params_ml(const vector<double> &vals_a,
				  const vector<double> &vals_b,
				  const vector<double> &probs) {
  const size_t lim = vals_a.size();
  if (workspace_vals_a.size() < lim) {
    workspace_vals_a.resize(lim);
    workspace_vals_b.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals_a[i] = log(vals_a[i]) + log(probs[i]);
    workspace_vals_b[i] = log(vals_b[i]) + log(probs[i]);
  }
  const double prob_sum = exp(split_distro_log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(split_distro_log_sum_log_vec(workspace_vals_a, lim))/prob_sum;
  params.back() = exp(split_distro_log_sum_log_vec(workspace_vals_b, lim))/prob_sum;
}

void
SkellamDistro::estimate_params_ml(const vector<double> &vals_a,
                                  const vector<double> &vals_b,
                                  const vector<double> &scales,
                                  const vector<double> &probs) {

  const size_t lim = vals_a.size();
  if (workspace_vals_a.size() < lim) {
    workspace_vals_a.resize(lim);
    workspace_vals_b.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i)
  {
      workspace_probs[i] = log(probs[i]) + log(scales[i]);
      workspace_vals_a[i] = log(vals_a[i]) + log(probs[i]);
      workspace_vals_b[i] = log(vals_b[i]) + log(probs[i]);
  }
  const double prob_sum = exp(split_distro_log_sum_log_vec(workspace_probs, lim));
  params.front() = exp(split_distro_log_sum_log_vec(workspace_vals_a, lim))/prob_sum;
  params.back() = exp(split_distro_log_sum_log_vec(workspace_vals_b, lim))/prob_sum;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

double
GaussianDistro::sample() const {
  assert(params.size() == 2);
  return gsl_ran_gaussian(SplitDistro_::rng, params.back()) + params.front();
}

double 
GaussianDistro::log_likelihood(const double val) const {
//   if (std::fabs(val) > 120)
//     return -40;
//   const double r = -(params.front() + params.back()) + 
//     (val/2.0)*(log(params.front()) - log(params.back())) +
//     log(gsl_sf_bessel_In(static_cast<int>(val), 
// 			 2*sqrt(params.front()*params.back())));
  if (gsl_ran_gaussian_pdf(val - params.front(), params.back()) > 1e-30) 
    return log(gsl_ran_gaussian_pdf(val - params.front(), params.back()));
  else
    return -30.0;
}

double 
GaussianDistro::log_likelihood(const double val,
                                   const double scale) const {

    const double mean = params[0] * scale;
    const double var = params[1] * pow(scale, 2.0);
    const double p = gsl_ran_gaussian_pdf(val - mean, var);
    return (p > 1e-30) ? p : -30.0;
}

GaussianDistro::GaussianDistro(const GaussianDistro &rhs) : 
  SplitDistro_(rhs.params) {}

GaussianDistro&
GaussianDistro::operator=(const GaussianDistro &rhs) {
  if (this != &rhs) {
    this->SplitDistro_::params = rhs.SplitDistro_::params;
    rng = gsl_rng_alloc(gsl_rng_default);
  }
  return *this;
}

void
GaussianDistro::estimate_params_ml(const vector<double> &vals_a, 
				  const vector<double> &vals_b) {
  assert(vals_a.size() == vals_b.size());
  const size_t lim = vals_a.size();

  vector<double> vals_c(vals_a);
  for (size_t i = 0; i < vals_c.size(); ++i)
    vals_c[i] -= vals_b[i];
  params.front() = accumulate(vals_c.begin(), vals_c.end(), 0.0)/lim;
  params.back() = gsl_stats_sd_m(&vals_c.front(), 1, vals_c.size(), params.front());
//  cerr << params.front() << "\t" << params.back() << endl;
}

void
GaussianDistro::estimate_params_ml(const vector<double> &vals_a,
				  const vector<double> &vals_b,
				  const vector<double> &probs) {
  const size_t lim = vals_a.size();
  if (workspace_vals_a.size() < lim) {
    workspace_vals_a.resize(lim);
    workspace_vals_b.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals_a[i] = log(vals_a[i]) + log(probs[i]);
    workspace_vals_b[i] = log(vals_b[i]) + log(probs[i]);
  }

  vector<double> vals_c(vals_a);
  for (size_t i = 0; i < vals_c.size(); ++i)
    vals_c[i] -= vals_b[i];
  params.front() = gsl_stats_wmean(&probs.front(), 1,
				   &vals_c.front(), 1, lim);
  params.back() = gsl_stats_wsd_m(&probs.front(), 1,
				  &vals_c.front(), 1, lim, params.front());
//  cerr << params.front() << "\t" << params.back() << endl;
}

void
GaussianDistro::estimate_params_ml(const vector<double> &vals_a,
                                   const vector<double> &vals_b,
                                   const vector<double> &scales,
                                   const vector<double> &probs) {

    cerr << "GaussianDistro::estimate_params_ml:  "
         << "Using an untested function" << endl;

  const size_t lim = vals_a.size();
  if (workspace_vals_a.size() < lim) {
    workspace_vals_a.resize(lim);
    workspace_vals_b.resize(lim);
    workspace_probs.resize(lim);
  }
  for (size_t i = 0; i < lim; ++i) {
      workspace_probs[i] = log(probs[i]);
    workspace_vals_a[i] = log(vals_a[i] / scales[i]) + log(probs[i]);
    workspace_vals_b[i] = log(vals_b[i] / scales[i]) + log(probs[i]);
  }

  vector<double> vals_c(vals_a);
  for (size_t i = 0; i < vals_c.size(); ++i)
      vals_c[i] = (vals_a[i] - vals_b[i]) / scales[i];
  params.front() = gsl_stats_wmean(&probs.front(), 1,
				   &vals_c.front(), 1, lim);
  params.back() = gsl_stats_wsd_m(&probs.front(), 1,
				  &vals_c.front(), 1, lim, params.front());
//  cerr << params.front() << "\t" << params.back() << endl;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void
NegBinomDiffDistro::set_helpers() {
  r1 = 1/params[1];
  p1 = r1/(r1 + params[0]);
  
  r2 = 1/params[3];
  p2 = r2/(r2 + params[2]);
  
  q1 = 1 - p1;
  q2 = 1 - p2;
  c = log(q1) + log(q2);
  
  r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2 = 
    (r1*log(p1) - gsl_sf_lngamma(r1) + r2*log(p2) - gsl_sf_lngamma(r2));
  //TODO: should check that these are valid!!!
}

double
NegBinomDiffDistro::sample() const {
  assert(params.size() == 4);
  const double mu1 = params[0];
  const double one_over_alpha1 = 1/params[1];
  const double mu2 = params[2];
  const double one_over_alpha2 = 1/params[3];
  return (int(gsl_ran_negative_binomial(SplitDistro_::rng, 
					one_over_alpha1/(one_over_alpha1 + mu1), 
					one_over_alpha1)) - 
	  int(gsl_ran_negative_binomial(SplitDistro_::rng, 
					one_over_alpha2/(one_over_alpha2 + mu2), 
					one_over_alpha2)));
}

const double NegBinomDiffDistro::tolerance = 1e-160;

static inline double
log_sum_log(const double p, const double q) {
  // Never pass a '0' into this function!!!!
  return ((p > q) ? p + log(1.0 + exp(q - p)) : q + log(1.0 + exp(p - q)));
}

inline static double
log_nbd_pdf(const unsigned int k, const double p, double n) {
  const double f = gsl_sf_lngamma(k + n);
  const double a = gsl_sf_lngamma(n);
  const double b = gsl_sf_lngamma(k + 1.0);
  return (f-a-b) + n*log(p) + k*log(1 - p);
}

double 
NegBinomDiffDistro::log_likelihood(const double k) const {
  const int k_int = static_cast<int>(k);
  if (ll_hash.find(k_int) == ll_hash.end()) {
    
    double sum = 0;
    double prev = -1;
    double n = std::max(0.0, -k);
    
    while (fabs((prev - sum)/sum) > tolerance) {
      prev = sum;
      const double log_form = 
 	(gsl_sf_lngamma(k + n + r1) - 
 	 gsl_sf_lnfact(static_cast<size_t>(k + n)) + 
 	 gsl_sf_lngamma(n + r2) - 
 	 gsl_sf_lnfact(static_cast<size_t>(n))) + n*c;
      sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
      ++n;
    }
    sum += k*log(q1) + r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2;
    ll_hash[k_int] = sum;
  }
  return ll_hash[k_int];
}

double 
NegBinomDiffDistro::log_likelihood(const double k,
                                   const double scale) const {

  const double scaled_p1 = r1 / (r1 + params[0] * scale);
  const double scaled_p2 = r2 / (r2 + params[2] * scale);
  const double scaled_log_q1_log_q2 = log(1 - scaled_p1)
      + log(1 - scaled_p2);
  const double scaled_r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2 =
      (r1 * log(scaled_p1) - gsl_sf_lngamma(r1)
       + r2 * log(scaled_p2) - gsl_sf_lngamma(r2));

  double sum = 0;
  double prev = -1;
  double n = std::max(0.0, -k);
    
  while (fabs((prev - sum)/sum) > tolerance)
  {
      prev = sum;
      const double log_form = 
          (gsl_sf_lngamma(k + n + r1) - 
           gsl_sf_lnfact(static_cast<size_t>(k + n)) + 
           gsl_sf_lngamma(n + r2) - 
           gsl_sf_lnfact(static_cast<size_t>(n))) + n * scaled_log_q1_log_q2;
      sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
      ++n;
  }
      
  sum += k * log(1 - scaled_p1) + scaled_r1_log_p1_M_lnG_r1_P_r2_log_p2_M_lnG_r2;
  return sum;
}


NegBinomDiffDistro::NegBinomDiffDistro(const NegBinomDiffDistro &rhs) : 
  SplitDistro_(rhs.params) {
  set_helpers();
}

NegBinomDiffDistro&
NegBinomDiffDistro::operator=(const NegBinomDiffDistro &rhs) {
  if (this != &rhs) {
    this->SplitDistro_::params = rhs.SplitDistro_::params;
    rng = gsl_rng_alloc(gsl_rng_default);
    set_helpers();
  }
  return *this;
}


static inline double
score_first_term(const vector<double> &v_hist, const double mu, const double alpha) {
  double sum = 0;
  const size_t lim = v_hist.size();
  for (size_t i = 0; i < lim; ++i)
    if (v_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
	inner_sum += j/(1.0 + alpha*j);
      sum += v_hist[i]*inner_sum;
    }
  return sum;
}

static inline double
alpha_score_function(const vector<double> &vals_hist, const double mu, 
		     const double alpha, const double vals_count) {
  const double one_plus_alpha_mu = 1 + alpha*mu;
  return (score_first_term(vals_hist, mu, alpha)/vals_count + 
	  (log(one_plus_alpha_mu)/alpha - mu)/alpha);
}

static const double max_allowed_alpha = 10;
static const double min_allowed_alpha = 1e-5;
static const double alpha_allowed_error = 1e-10;

static void
estimate_params_pair_ml(const vector<double> &vals, double &mu, double &alpha) {
  const size_t lim = vals.size();
  // This is the mu
  mu = std::accumulate(vals.begin(), vals.end(), 0.0)/lim;
  
  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.end());
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < lim; ++i)
    ++vals_hist[static_cast<size_t>(vals[i])];
  
  const double vals_count = lim;
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0) 
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;
  // METHOD OF MOMENTS:
  //   const size_t lim = vals.size();
  //   mu = std::accumulate(vals.begin(), vals.begin() + lim, 0.0)/lim;
  //   const double var = gsl_stats_variance_m(&vals.front(), 1, vals.size(), mu);
  //   const double r = (mu*mu)/(var - mu);
  //   alpha = max(min_alpha, 1/r);
}

static void
estimate_params_pair_ml(const vector<double> &vals, const vector<double> &probs,
			vector<double> &workspace_vals, vector<double> &workspace_probs,
			double &mu, double &alpha) {
  const size_t lim = vals.size();
  if (workspace_vals.size() < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }

  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }
  
  const double vals_count = exp(SplitDistro::log_sum_log_vec(workspace_probs, lim));
  mu = exp(SplitDistro::log_sum_log_vec(workspace_vals, lim))/vals_count;
  
  const double max_value = *std::max_element(vals.begin(), vals.end());
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < lim; ++i)
    vals_hist[static_cast<size_t>(vals[i])] += probs[i];
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;
  // METHOD OF MOMENTS:
  //   mu = exp(split_distro_log_sum_log_vec(workspace_vals, lim))/vals_count;
  //   const double var = gsl_stats_wvariance_m(&vals.front(), 1, &probs.front(), 1, 
  // 					   vals.size(), mu);
  //   const double r = (mu*mu)/(var - mu);
  //   alpha = max(min_alpha, 1/r);
}

void
NegBinomDiffDistro::estimate_params_ml(const vector<double> &vals_a, 
                                       const vector<double> &vals_b) {
    hq_estimate_params_ml(vals_a, vals_b);
}

void
NegBinomDiffDistro::estimate_params_ml(const vector<double> &vals_a,
				       const vector<double> &vals_b,
				       const vector<double> &probs) {
    hq_estimate_params_ml(vals_a, vals_b, probs);
}


////// andrew's orignial version
void
NegBinomDiffDistro::andrew_estimate_params_ml(const vector<double> &vals_a, 
                                       const vector<double> &vals_b) {
  estimate_params_pair_ml(vals_a, params[0], params[1]);
  estimate_params_pair_ml(vals_b, params[2], params[3]);
  
  ll_hash.clear();
  set_helpers();
}


void
NegBinomDiffDistro::andrew_estimate_params_ml(const vector<double> &vals_a,
				       const vector<double> &vals_b,
				       const vector<double> &probs) {
  estimate_params_pair_ml(vals_a, probs, workspace_vals_a, workspace_probs,
			  params[0], params[1]);
  estimate_params_pair_ml(vals_b, probs, workspace_vals_b, workspace_probs,
			  params[2], params[3]);
  ll_hash.clear();
  set_helpers();
}
///////

/////// experimental code for gradient descent method 
inline void
my_dummy_func(const double &r)
{
    // used to make sure a tempory variable are not optimized out by compiler 
    return;
}

double
local_log_likelihood(const vector<double> &vals, const vector<double> &probs,
                     const double &mu1, const double &alpha1,
                     const double &mu2, const double &alpha2)
{
    static const double tolerance = 1e-160;
    const double r1 = 1 / alpha1;
    const double p1 = r1 / (r1 + mu1);
    const double log_q1 = log(1 - p1);

    const double r2 = 1 / alpha2;
    const double p2 = r2 / (r2 + mu2);
    
    const double fixed_term = (r1*log(p1) - gsl_sf_lngamma(r1)
                               + r2*log(p2) - gsl_sf_lngamma(r2));

    const double log_q1_log_q2 = log(1 - p1) + log(1 - p2);

    unordered_map<int, double> ll_hash;
    
    double llh = 0;
    for (size_t i = 0; i < vals.size(); ++i)
    {
        const int val_int = static_cast<int>(vals[i]);
        if (ll_hash.find(val_int) == ll_hash.end())
        {
            double sum = 0;
            double prev = -1;
            double n = std::max(0.0, - vals[i]);
            
            while (fabs((prev - sum)/sum) > tolerance) {
                prev = sum;
                const double log_form = 
                    (gsl_sf_lngamma(val_int + n + r1) - 
                     gsl_sf_lnfact(static_cast<size_t>(val_int + n)) + 
                     gsl_sf_lngamma(n + r2) - 
                     gsl_sf_lnfact(static_cast<size_t>(n))) + n * log_q1_log_q2;
                sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
                ++n;
            }
            sum += val_int * log_q1 + fixed_term;
            ll_hash[val_int] = sum;
        }
        llh += probs[i] * ll_hash[val_int];
    }

    return llh;
}

double
local_log_likelihood(const vector<double> &vals,
                     const double &mu1, const double &alpha1,
                     const double &mu2, const double &alpha2)
{
    const vector<double> probs(vals.size(), 1);
    return local_log_likelihood(vals, probs, mu1, alpha1,
                                mu2, alpha2);
}


void
NegBinomDiffDistro::hq_estimate_params_ml(const vector<double> &vals_a, 
                                          const vector<double> &vals_b)
{
    const vector<double> probs(vals_a.size(), 1);
    hq_estimate_params_ml(vals_a, vals_b, probs);
}

void 
NegBinomDiffDistro::hq_estimate_params_ml(const std::vector<double> &vals_a,
                                          const std::vector<double> &vals_b,
                                          const std::vector<double> &probs)
{

//  static const double max_allow_error = 1e-10;
    static const double max_allow_llh_diff = 1e-12;
//    static const double max_delta_norm = 1e-20;
    static const double sqrt_error_machine = 1e-10;

    double global_scale = 1.0 / 16;
    
    vector<double> vals(vals_a.size());
    for (size_t i = 0; i < vals.size(); ++i)
        vals[i] = vals_a[i] - vals_b[i];
    
    double mu1, alpha1, mu2, alpha2;

    // get initial value
    estimate_params_pair_ml(vals_a, probs,
                            workspace_vals_a, workspace_probs,
                            mu1, alpha1);
    estimate_params_pair_ml(vals_b, probs,
                            workspace_vals_b, workspace_probs,
                            mu2, alpha2);
    
    // newton gradient
    double temp(0);
    double h(0);
    
    double mu1_grad(0), alpha1_grad(0), mu2_grad(0), alpha2_grad(0);
    double prev_llh = local_log_likelihood(vals, probs, mu1, alpha1, mu2, alpha2);
    double llh_diff = std::numeric_limits<double>::max();
    double delta_norm = std::numeric_limits<double>::max();
    double curr_llh = 0;
    
    size_t iter(0);
    while (fabs(llh_diff / curr_llh) > max_allow_llh_diff)
    {
        ++iter;
        mu1 += mu1_grad;
        alpha1 += alpha1_grad;
        mu2 += mu2_grad;
        alpha2 += alpha2_grad;
        
        // compute gradient
        h = sqrt_error_machine * mu1; // stepwise to be determined wisely
        temp = h + mu1;
        my_dummy_func(temp);
        h = temp - mu1;
        mu1_grad  = (local_log_likelihood(vals, probs, mu1 + h, alpha1, mu2, alpha2)
                     -  local_log_likelihood(vals, probs, mu1 - h, alpha1, mu2, alpha2))
            / h / 2;
        
        h = sqrt_error_machine * alpha1; // stepwise to be determined wisely
        temp = h + alpha1;
        my_dummy_func(temp);
        h = temp - alpha1;
        alpha1_grad  = (local_log_likelihood(vals, probs, mu1, alpha1 + h, mu2, alpha2)
                        - local_log_likelihood(vals, probs, mu1, alpha1 - h, mu2, alpha2))
            / h / 2;
        
        h = sqrt_error_machine * mu2; // stepwise to be determined wisely
        temp = h + mu2;
        my_dummy_func(temp);
        h = temp - mu2;
        mu2_grad  = (local_log_likelihood(vals, probs, mu1, alpha1, mu2 + h, alpha2)
                     - local_log_likelihood(vals, probs, mu1, alpha1, mu2 - h, alpha2))
            / h / 2;
        
        h = sqrt_error_machine * alpha2; // stepwise to be determined wisely
        temp = h + alpha2;
        my_dummy_func(temp);
        h = temp - alpha2;
        alpha2_grad = (local_log_likelihood(vals, probs, mu1, alpha1, mu2, alpha2 + h)
                       - local_log_likelihood(vals, probs, mu1, alpha1, mu2, alpha2 - h))
            / h / 2;
        
        double local_scale = min(min(mu1 / fabs(mu1_grad + 1e-20),
                                     alpha1 / fabs(alpha1_grad + 1e-20)),
                                 min(mu2 / fabs(mu2_grad + 1e-20),
                                     alpha2 / fabs(alpha2_grad + 1e-20)));
        mu1_grad *= local_scale * global_scale;
        alpha1_grad *= local_scale * global_scale;
        mu2_grad *= local_scale * global_scale;
        alpha2_grad *= local_scale * global_scale;
    
        curr_llh = local_log_likelihood(vals, probs,
                                        mu1 + mu1_grad, alpha1 + alpha1_grad,
                                        mu2 + mu2_grad, alpha2 + alpha2_grad);
        llh_diff = curr_llh - prev_llh;
        
        while (llh_diff < - 1e-100) // if the move can not lead to greater likelihood
            // cancel this remove and try a smaller step
        {
            mu1_grad /= 2;
            alpha1_grad /= 2;
            mu2_grad /= 2;
            alpha2_grad /= 2;
            
            curr_llh = local_log_likelihood(vals, probs,
                                            mu1 + mu1_grad, alpha1 + alpha1_grad,
                                            mu2 + mu2_grad, alpha2 + alpha2_grad);
            llh_diff = curr_llh - prev_llh;
        }
        prev_llh = curr_llh;
        delta_norm = pow(mu1_grad, 2.0) + pow(alpha1_grad, 2.0) +
            pow(mu2_grad, 2.0) + pow(alpha2_grad, 2.0);
    }

    params[0] = mu1;
    params[1] = alpha1;
    params[2] = mu2;
    params[3] = alpha2;
    
    ll_hash.clear();
    set_helpers();
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

double
llh_deriative_rt_alpha(const vector<double> &vals,
                       const vector<double> &scales,
                       const vector<double> &probs,
                       const vector<double> &vals_hist,
                       const double mu,
                       const double alpha);
// defined in Distro.cpp
// {
//     const double first_term = score_fun_first_term(vals_hist, mu, alpha);
    
//     const double mu_times_alpha = mu * alpha;
//     const double alpha_inverse = 1 / alpha;
//     const double alpha_square_inverse = pow(alpha_inverse, 2.0);
    
//     double second_term = 0;
//     for (size_t i = 0; i < vals.size(); ++i)
//     {
//         const double one_plus_extra = 1 + scales[i] * mu_times_alpha;
//         second_term += probs[i] *
//             (alpha_square_inverse * log(one_plus_extra) -
//              scales[i] * mu * (alpha_inverse + vals[i]) / one_plus_extra);
//     }

//     return first_term + second_term;
// }

void 
estimate_params_pair_ml(const std::vector<double> &vals,
                        const std::vector<double> &scales,
                        const std::vector<double> &probs,
                        vector<double> &workspace_vals,
                        vector<double> &workspace_probs,
                        double &mu, double &alpha)
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
    mu = exp(split_distro_log_sum_log_vec(workspace_vals, lim)) /
        exp(split_distro_log_sum_log_vec(workspace_probs, lim));

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
    alpha = a_mid;
}

double
local_log_likelihood(const vector<double> &vals,
                     const vector<double> &scales,
                     const vector<double> &probs,
                     const double &mu1, const double &alpha1,
                     const double &mu2, const double &alpha2)
{
    static const double tolerance = 1e-160;

    double llh = 0;
    for (size_t i = 0; i < vals.size(); ++i)
    {
        const double r1 = 1 / alpha1;
        const double p1 = r1 / (r1 + mu1 * scales[i]);
        const double log_q1 = log(1 - p1);
        
        const double r2 = 1 / alpha2;
        const double p2 = r2 / (r2 + mu2 * scales[i]);
        
        const double fixed_term = (r1*log(p1) - gsl_sf_lngamma(r1)
                                   + r2*log(p2) - gsl_sf_lngamma(r2));
        
        const double log_q1_log_q2 = log(1 - p1) + log(1 - p2);

        const int val_int = static_cast<int>(vals[i]);

        double sum = 0;
        double prev = -1;
        double n = std::max(0.0, - vals[i]);
            
        while (fabs((prev - sum)/sum) > tolerance)
        {
            prev = sum;
            const double log_form = 
                (gsl_sf_lngamma(val_int + n + r1) - 
                 gsl_sf_lnfact(static_cast<size_t>(val_int + n)) + 
                 gsl_sf_lngamma(n + r2) - 
                 gsl_sf_lnfact(static_cast<size_t>(n))) + n * log_q1_log_q2;
            sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
                ++n;
        }
        sum += val_int * log_q1 + fixed_term;

        llh += probs[i] * sum;
    }

    return llh;
}


void 
NegBinomDiffDistro::estimate_params_ml(const std::vector<double> &vals_a,
                                       const std::vector<double> &vals_b,
                                       const std::vector<double> &scales,
                                       const std::vector<double> &probs)
{
    estimate_params_pair_ml(vals_a, scales, probs,
                            workspace_vals_a, workspace_probs,
                            params[0], params[1]);
    estimate_params_pair_ml(vals_b, scales, probs,
                            workspace_vals_b, workspace_probs,
                            params[2], params[3]);

//     vector<double> vas(vals_a), vbs(vals_b);
    
//     for (size_t i = 0; i < vas.size(); ++i)
//     {
//         double c = 0;
        
//         vas[i] /= scales[i];
//         c = floor(vas[i]);
//         vas[i] = c + static_cast<int>(vas[i] - c > 0.5);

//         vbs[i] /= scales[i];
//         c = floor(vbs[i]);
//         vbs[i] = c + static_cast<int>(vbs[i] - c > 0.5);
//     }

//     estimate_params_ml(vas, vbs, probs);
     
//////// the exact implementation of MLE methods with scales data
//     however this code is very slow
//     static const double max_allow_llh_diff = 1e-8;
//     static const double max_delta_norm = 1e-20;
//     static const double sqrt_error_machine = 1e-10;

//     double global_scale = 1.0 / 16;
//     vector<double> vals(vals_a.size());
//     for (size_t i = 0; i < vals.size(); ++i)
//         vals[i] = vals_a[i] - vals_b[i];
    
//     double mu1, alpha1, mu2, alpha2;

//     // get initial value
//     estimate_params_pair_ml(vals_a, scales, probs,
//                             workspace_vals_a, workspace_probs,
//                             mu1, alpha1);
//     estimate_params_pair_ml(vals_b, scales, probs,
//                             workspace_vals_b, workspace_probs,
//                             mu2, alpha2);
    
//     // newton gradient
//     double temp(0);
//     double h(0);
    
//     double mu1_grad(0), alpha1_grad(0), mu2_grad(0), alpha2_grad(0);
//     double prev_llh = local_log_likelihood(vals, scales, probs, mu1, alpha1, mu2, alpha2);
//     double llh_diff = std::numeric_limits<double>::max();
//     double delta_norm = std::numeric_limits<double>::max();
//     double curr_llh = 0;
    
//     size_t iter(0);
//     while (fabs(llh_diff / curr_llh) > max_allow_llh_diff &&
//            delta_norm > max_delta_norm)
//     {
//         ++iter;
//         mu1 += mu1_grad;
//         alpha1 += alpha1_grad;
//         mu2 += mu2_grad;
//         alpha2 += alpha2_grad;
        
//         // compute gradient
//         h = sqrt_error_machine * mu1; // stepwise to be determined wisely
//         temp = h + mu1;
//         my_dummy_func(temp);
//         h = temp - mu1;
//         mu1_grad  = (local_log_likelihood(vals, scales, probs, mu1 + h, alpha1, mu2, alpha2)
//                      -  prev_llh) / h;
        
//         h = sqrt_error_machine * alpha1; // stepwise to be determined wisely
//         temp = h + alpha1;
//         my_dummy_func(temp);
//         h = temp - alpha1;
//         alpha1_grad  = (local_log_likelihood(vals, scales, probs, mu1, alpha1 + h, mu2, alpha2)
//                         - prev_llh) / h;
        
//         h = sqrt_error_machine * mu2; // stepwise to be determined wisely
//         temp = h + mu2;
//         my_dummy_func(temp);
//         h = temp - mu2;
//         mu2_grad  = (local_log_likelihood(vals, scales, probs, mu1, alpha1, mu2 + h, alpha2)
//                      - prev_llh) / h;
        
//         h = sqrt_error_machine * alpha2; // stepwise to be determined wisely
//         temp = h + alpha2;
//         my_dummy_func(temp);
//         h = temp - alpha2;
//         alpha2_grad = (local_log_likelihood(vals, scales, probs, mu1, alpha1, mu2, alpha2 + h)
//                        - prev_llh) / h;
        
//         double local_scale = min(min(mu1 / fabs(mu1_grad + 1e-20),
//                                      alpha1 / fabs(alpha1_grad + 1e-20)),
//                                  min(mu2 / fabs(mu2_grad + 1e-20),
//                                      alpha2 / fabs(alpha2_grad + 1e-20)));
//         mu1_grad *= local_scale * global_scale;
//         alpha1_grad *= local_scale * global_scale;
//         mu2_grad *= local_scale * global_scale;
//         alpha2_grad *= local_scale * global_scale;
        
//         curr_llh = local_log_likelihood(vals, scales, probs,
//                                         mu1 + mu1_grad, alpha1 + alpha1_grad,
//                                         mu2 + mu2_grad, alpha2 + alpha2_grad);
//         llh_diff = curr_llh - prev_llh;
        
//         while (llh_diff < - 1e-100) // if the move can not lead to greater likelihood
//             // cancel this remove and try a smaller step
//         {
//             mu1_grad /= 2;
//             alpha1_grad /= 2;
//             mu2_grad /= 2;
//             alpha2_grad /= 2;
            
//             curr_llh = local_log_likelihood(vals, scales, probs,
//                                             mu1 + mu1_grad, alpha1 + alpha1_grad,
//                                             mu2 + mu2_grad, alpha2 + alpha2_grad);
//             llh_diff = curr_llh - prev_llh;
//         }
//         prev_llh = curr_llh;
//         delta_norm = pow(mu1_grad, 2.0) + pow(alpha1_grad, 2.0) +
//             pow(mu2_grad, 2.0) + pow(alpha2_grad, 2.0);
//     }

//     params[0] = mu1;
//     params[1] = alpha1;
//     params[2] = mu2;
//     params[3] = alpha2;
    
//     ll_hash.clear();
//     set_helpers();
}


///////


// //////////  experimental code for BFGS method: hopefully faster
// inline void
// colvec_rowvec(const vector<double> &colvec,
//               const vector<double> &rowvec,
//               vector< vector<double> > &m)
// {
//     for (size_t i = 0; i < colvec.size(); ++i)
//         for (size_t j = 0; j < rowvec.size(); ++j)
//             m[i][j] = colvec[i] * rowvec[j];
// }

// inline double
// rowvec_colvec(const vector<double> &rowvec,
//               const vector<double> &colvec)
// {
//     return std::inner_product(rowvec.begin(), rowvec.end(),
//                               colvec.begin(), 0.0);
// }

// inline void
// matrix_colvec(const vector< vector<double> > &m,
//               const vector<double> &colvec,
//               vector<double> &re_colvec)
// {
//     for (size_t i = 0; i < m.size(); ++i)
//         re_colvec[i] = rowvec_colvec(m[i], colvec);
// }

// inline void
// rowvec_matrix(const vector<double> &rowvec,
//               const vector< vector<double> > &m,
//               vector<double> &re_rowvec)
// {
//     std::fill(re_rowvec.begin(), re_rowvec.end(), 0.0);
    
//     const size_t n = rowvec.size();
//     for (size_t i = 0; i < n; ++i)
//         for (size_t j = 0; j < n; ++j)
//             re_rowvec[i] += rowvec[j] * m[j][i];
// }

// inline double
// rowvec_matrix_colvec(const vector<double> &rowvec,
//                      const vector< vector<double> > &m,
//                      const vector<double> &colvec)
// {
//     static vector<double> v(rowvec.size(), 0);
//     rowvec_matrix(rowvec, m, v);
//     return rowvec_colvec(v, colvec);
// }

// void 
// NegBinomDiffDistro::bfgs_estimate_params_ml(const std::vector<double> &vals_a, 
//                                             const std::vector<double> &vals_b)
// {
//     const vector<double> probs(vals_a.size(), 1);
//     bfgs_estimate_params_ml(vals_a, vals_b, probs);
// }

// void 
// NegBinomDiffDistro::bfgs_estimate_params_ml(const std::vector<double> &vals_a, 
//                                             const std::vector<double> &vals_b,
//                                             const vector<double> &probs)
// {
//     static const double max_allow_llh_diff = 1e-10;
//     static const double sqrt_error_machine = 1e-8;
//     const size_t params_n = 4;
    
//     double global_scale = 1.0 / 8;
//     double local_scale = 1;

//     vector<double> vals(vals_a.size());
//     for (size_t i = 0; i < vals.size(); ++i)
//         vals[i] = vals_a[i] - vals_b[i];

//     // variables and inilization
//     vector<double> xs(params_n, 0); // parameters
//     vector<double> delta_xs(params_n, 0);
//     vector<double> gradients(params_n, 0);
//     vector<double> prev_gradients(params_n, 0);
//     vector<double> delta_grads(params_n, 0);
//     vector< vector<double> > hi(params_n, vector<double>(params_n, 0)); // Hemesian inverse

//     estimate_params_pair_ml(vals_a, probs,
//                             workspace_vals_a, workspace_probs,
//                             xs[0], xs[1]);
//     estimate_params_pair_ml(vals_b, probs,
//                             workspace_vals_b, workspace_probs,
//                             xs[2], xs[3]);
//     for (size_t i = 0; i < params_n; ++i) hi[i][i] = 1; // start from identity matrix
    
//     // temporary var
//     double coeff1(0), coeff2(0);
//     vector<double> tmp_v1(params_n, 0), tmp_v2(params_n, 0);
//     double h(0), temp(0);
    

//     // get initial values ready for big while loop
//     double prev_llh =  - local_log_likelihood(vals, probs, xs[0], xs[1], xs[2], xs[3]);
//     double llh_diff = std::numeric_limits<double>::max();
//     double curr_llh = 0;

//     for (size_t i = 0; i < params_n; ++i) // initial gradient
//     {
//         h = sqrt_error_machine * xs[i]; // stepwise to be determined wisely
//         temp = h + xs[i]; // avoid round off error
//         my_dummy_func(temp);
//         h = temp - xs[i];
//         xs[i] += h;
//         gradients[i]  = (- local_log_likelihood(vals, probs, xs[0], xs[1], xs[2], xs[3])
//                          + prev_llh) / h;
//         xs[i] -= h;
//     }
    
//     while (fabs(llh_diff / curr_llh) > max_allow_llh_diff && true)
//     {
//         // determine direction and step size
//         matrix_colvec(hi, gradients, delta_xs);
//         local_scale = numeric_limits<double>::max();
//         for (size_t i = 0; i < params_n; ++i)
//         {
//             const double alpha = fabs(xs[i] / (delta_xs[i] + 1e-20));
//             local_scale = (local_scale < alpha) ? local_scale : alpha;
//         }
//         for (size_t i = 0; i < params_n; ++i)
//             delta_xs[i] *= -local_scale * global_scale;
        
//         // determine whether step is appropriviate
//         curr_llh = - local_log_likelihood(vals, probs,
//                                           xs[0] + delta_xs[0], xs[1] + delta_xs[1],
//                                           xs[2] + delta_xs[2], xs[3] + delta_xs[3]);
//         llh_diff = curr_llh - prev_llh;
        
//         // if the move can not lead to greater likelihood
//         // cancel this remove and try a smaller step
//         while (llh_diff >  0)
//         {
//             global_scale /= 2;
//             for (size_t i = 0; i < params_n; ++i)
//                 delta_xs[i] /= 2;
            
//             curr_llh = - local_log_likelihood(vals, probs,
//                                               xs[0] + delta_xs[0], xs[1] + delta_xs[1],
//                                               xs[2] + delta_xs[2], xs[3] + delta_xs[3]);
//             llh_diff = curr_llh - prev_llh;
//         }

// /* debug code  */
//         cerr << "check:" << endl
//              << xs[0] << "\t" << xs[1] << "\t"<< xs[2] << "\t"<< xs[3] << endl
//              << delta_xs[0] << "\t" << delta_xs[1] << "\t"<< delta_xs[2] << "\t"<< delta_xs[3] << endl
//              << -prev_llh << "\t" << curr_llh << endl;
// // */        

//         prev_llh = curr_llh;


//         // update estimated values
//         for (size_t i  = 0; i < params_n; ++i)
//             xs[i] += delta_xs[i];
        
//         // update gradients
//         std::copy(gradients.begin(), gradients.end(), prev_gradients.begin());

//         for (size_t i = 0; i < params_n; ++i)
//         {
//             h = sqrt_error_machine * xs[i]; // stepwise to be determined wisely
//             temp = h + xs[i]; // avoid round off error
//             my_dummy_func(temp);
//             h = temp - xs[i];
//             xs[i] += h;
//             gradients[i]  = (local_log_likelihood(vals, probs,
//                                                   xs[0], xs[1], xs[2], xs[3])
//                              - prev_llh) / h;
//             xs[i] -= h;
//         }
        
//         for (size_t i = 0; i < params_n; ++i)
//             delta_grads[i] = gradients[i] - prev_gradients[i];

//         // update Hessian Inverse
//         coeff2 = 1 / rowvec_colvec(delta_xs, delta_grads);
//         coeff1 = 0 +
//             rowvec_matrix_colvec(delta_grads, hi, delta_grads) * pow(coeff2, 2.0);
//         matrix_colvec(hi, delta_grads, tmp_v1);
//         rowvec_matrix(delta_grads, hi, tmp_v2);
        
//         for (size_t i = 0; i < params_n; ++i)
//             for (size_t j = 0; j < params_n; ++j)
//                 hi[i][j] += delta_xs[i] * delta_xs[j] * coeff1
//                     - (tmp_v1[i]*delta_xs[j] + delta_xs[i]*tmp_v2[j]) * coeff2;
//     }
    
//     // store result
//     std::copy(xs.begin(), xs.end(), params.begin());
//     ll_hash.clear();
//     set_helpers();
// }

// //////////


// /////////// call gsl version of BFGS minimizer
// double 
// minus_llh(const gsl_vector *xs, void *ps)
// {
//     const double *arr = (double *)ps;
//     const int n = (int)(arr[0]);
//     const double *vals = arr + 1;
//     const double *probs = arr + n;
    
//     const double mu1 = gsl_vector_get(xs, 0);
//     const double alpha1 = gsl_vector_get(xs, 1) > 1e-5 ? gsl_vector_get(xs, 1) : 1e-5;
//     const double mu2 = gsl_vector_get(xs, 2);
//     const double alpha2 = gsl_vector_get(xs, 2) > 1e-5 ? gsl_vector_get(xs, 2) : 1e-5;
    
//     const double tolerance = 1e-160;
//     const double r1 = 1 / alpha1;
//     const double p1 = r1 / (r1 + mu1);
//     const double log_q1 = log(1 - p1);

//     const double r2 = 1 / alpha2;
//     const double p2 = r2 / (r2 + mu2);
    
//     const double fixed_term = (r1*log(p1) - gsl_sf_lngamma(r1)
//                                + r2*log(p2) - gsl_sf_lngamma(r2));

//     const double log_q1_log_q2 = log(1 - p1) + log(1 - p2);

//     double ll_table[4000];
//     const int offset = 2000;

//     int i;
//     for (i = 0; i < 4000; ++i) ll_table[i] = 10;
    
//     double llh = 0;
//     for (i = 0; i < n; ++i)
//     {
//         double val = *(vals + i);
//         double prob = *(probs + i);
//         int val_int = (int)(val);
//         if (ll_table[val_int + offset] > 0)
//         {
//             double sum = 0;
//             double prev = -1;
//             double ui = (0.0 > - val) ? 0 : -val;
            
//             while (fabs((prev - sum)/sum) > tolerance)
//             {
//                 prev = sum;
//                 double log_form = 
//                     (gsl_sf_lngamma(val_int + ui + r1) - 
//                      gsl_sf_lnfact((int)(val_int + ui)) + 
//                      gsl_sf_lngamma(ui + r2) - 
//                      gsl_sf_lnfact((int)(ui))) + ui * log_q1_log_q2;
//                 sum = ((sum == 0) ? log_form : log_sum_log(sum, log_form));
//                 ++ui;
//             }
//             sum += val_int * log_q1 + fixed_term;
//             ll_table[val_int + offset] = sum;
//         }
//         llh += prob * ll_table[val_int + offset];
//     }

// //  cerr << "llh done" << endl;

//     return -llh;
// }


// void
// minus_llh_df(const gsl_vector *xs, void *ps,
//              gsl_vector *df)
// {
//     const double sqrt_error_machine = 1e-6;
//     size_t i;
//     gsl_vector *tv1 = gsl_vector_alloc(4);
//     gsl_vector *tv2 = gsl_vector_alloc(4);
//     gsl_vector_memcpy(tv1, xs);
//     gsl_vector_memcpy(tv2, xs);
    
//     for (i = 0; i < 4; ++i) // initial gradient
//     {
//         double x = gsl_vector_get(xs, i);
//         double h = sqrt_error_machine * x; // stepwise to be determined wisely
//         double temp = h + x; // avoid round off error
//         my_dummy_func(temp);
//         h = temp - x;
//         *gsl_vector_ptr(tv1, i) += h;
//         *gsl_vector_ptr(tv2, i) -= h;
//         double d  = (minus_llh(tv1, ps) - minus_llh(tv2, ps)) / h / 2;
//         *gsl_vector_ptr(df, i) = d;
//         *gsl_vector_ptr(tv1, i) -= h;
//         *gsl_vector_ptr(tv2, i) += h;
//     }
//     gsl_vector_free(tv1);
//     gsl_vector_free(tv2);
// }

// void
// minus_llh_fdf(const gsl_vector *xs, void * ps,
//               double *f, gsl_vector *df)
// {
//     *f = minus_llh(xs, ps);
//     minus_llh_df(xs, ps, df);
// }

// void
// NegBinomDiffDistro::gsl_bfgs_estimate_params_ml(const std::vector<double> &vals_a, 
//                                                 const std::vector<double> &vals_b)
// {
//     const vector<double> probs(vals_a.size(), 1);
//     gsl_bfgs_estimate_params_ml(vals_a, vals_b, probs);
// }

// void 
// NegBinomDiffDistro::gsl_bfgs_estimate_params_ml(const std::vector<double> &vals_a, 
//                                                 const std::vector<double> &vals_b,
//                                                 const vector<double> &probs)
// {
//     cerr << "check:" << "initialize minimizer"  << endl;
    
//     const size_t sz = probs.size();
//     double ps[sz * 2 + 1];
//     ps[0] = sz;
//     for (size_t i = 0; i < sz; ++i)
//     {
//         ps[i + 1] = vals_a[i] - vals_b[i];
//         ps[i + sz + 1] = probs[i];
//     }

//     estimate_params_pair_ml(vals_a, probs, workspace_vals_a, workspace_probs,
//                             params[0], params[1]);
//     estimate_params_pair_ml(vals_b, probs, workspace_vals_b, workspace_probs,
//                             params[2], params[3]);

//     gsl_vector * xs = gsl_vector_alloc(4);
//     gsl_vector_set(xs, 0, params[0]);
//     gsl_vector_set(xs, 1, params[1]);
//     gsl_vector_set(xs, 2, params[2]);
//     gsl_vector_set(xs, 3, params[3]);

//     gsl_multimin_function_fdf my_func;
//     my_func.n = 4;  /* number of function components */
//     my_func.f = &minus_llh;
//     my_func.df = &minus_llh_df;
//     my_func.fdf = &minus_llh_fdf;
//     my_func.params = (void *)ps;

//     const gsl_multimin_fdfminimizer_type *T;
//     gsl_multimin_fdfminimizer *optimizer;
    
//     T = gsl_multimin_fdfminimizer_steepest_descent;
//     optimizer = gsl_multimin_fdfminimizer_alloc(T, 4);
    
//     int iter = 0, status;
//     gsl_multimin_fdfminimizer_set(optimizer, &my_func, xs, 0.01, 1e-10);

//     cerr << "llh = " <<  - minus_llh(xs, (void *)ps) << endl;
    
//     cerr << "check:" << "run minimizer"  << endl;

//     do
//     {
//         iter++;
//         status = gsl_multimin_fdfminimizer_iterate(optimizer);

//         cerr << "check 1.0"  << endl;
        
//         if (status)
//             break;
        
//         status = gsl_multimin_test_gradient(optimizer->gradient, 1e-10);

// //*/
//         cerr << "check : " << endl
//              << *gsl_vector_ptr(optimizer->x, 0) << "\t"
//              << *gsl_vector_ptr(optimizer->x, 1) << "\t"
//              << *gsl_vector_ptr(optimizer->x, 2) << "\t"
//              << *gsl_vector_ptr(optimizer->x, 3) << endl 
//              << *gsl_vector_ptr(optimizer->dx, 0) << "\t"
//              << *gsl_vector_ptr(optimizer->dx, 1) << "\t"
//              << *gsl_vector_ptr(optimizer->dx, 2) << "\t"
//              << *gsl_vector_ptr(optimizer->dx, 3) << endl
//              << *gsl_vector_ptr(optimizer->gradient, 0) << "\t"          
//              << *gsl_vector_ptr(optimizer->gradient, 1) << "\t"
//              << *gsl_vector_ptr(optimizer->gradient, 2) << "\t"
//              << *gsl_vector_ptr(optimizer->gradient, 3) << endl
//              << optimizer->f << endl;
// //*/
        
        
//         if (status == GSL_SUCCESS)
//             cerr <<  "Minimum found" << endl;
//     }
//     while (status == GSL_CONTINUE && iter < 100);
    
//     cerr << "check:" << "finish minimizer"  << endl;

    
//     for (size_t i = 0; i < 4; ++i)
//         params[i] = *gsl_vector_ptr(xs, i);
//     ll_hash.clear();
//     set_helpers();

    
//     gsl_multimin_fdfminimizer_free(optimizer);    
//     gsl_vector_free(xs);
// }
// ///////////




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////
////// DISCRETE EMPIRICAL SPLIT DISTRIBUTION
//////

const double DiscEmpSplitDistro::MIN_PROB = 1e-20;
const double DiscEmpSplitDistro::smoothing_parameter = 10.0;

size_t
DiscEmpSplitDistro::find_bin(const vector<double> &bins, 
			     const double shifted_val) {
  const size_t bin = upper_bound(bins.begin(), bins.end(), shifted_val) - 
    bins.begin() - 1;
  assert(bin < bins.size() && bin >= 0);
  return bin;
}


void
DiscEmpSplitDistro::make_cumulative(std::vector<double> &shifted_vals) {
  for (size_t i = 1; i < shifted_vals.size(); ++i)
    shifted_vals[i] += shifted_vals[i - 1];
}


void
DiscEmpSplitDistro::make_hist(const vector<double> &shifted_data, 
			      const size_t n_classes, 
			      const double max_val, const double min_val,
			      vector<double> &hist) {
  const size_t n_vals = shifted_data.size();
  
  hist.resize(n_classes, 0);
  for (size_t i = 0; i < n_vals; ++i) {
    assert(shifted_data[i] >= 0 && shifted_data[i] <= n_classes);
    hist[static_cast<size_t>(shifted_data[i])]++; // Already shifted
  }
  hist.back() = 0;
}


void
DiscEmpSplitDistro::make_weighted_hist(const vector<double> &shifted_data,
				       const vector<double> &weights, 
				       const size_t n_classes, 
				       const double max_val, 
				       const double min_val,
				       vector<double> &hist) {
  
  const size_t n_vals = shifted_data.size();
  // fill the counts
  hist.resize(n_classes, 0.0);
  for (size_t i = 0; i < n_vals; ++i) {
    assert(shifted_data[i] >= 0 && shifted_data[i] <= n_classes);
    hist[static_cast<size_t>(shifted_data[i])] += weights[i]; // Already shifted
  }
  hist.back() = 0;
}


double
DiscEmpSplitDistro::sample() const {
  assert(!params.empty());
  const size_t id = find_bin(cumulative, gsl_rng_uniform(rng));
  assert(id >= min_val && id <= max_val);
  return id + min_val; // Need to shift the bin id
}

double 
DiscEmpSplitDistro::log_likelihood(const double unshifted_val) const {
  // need to shift
  const size_t index = static_cast<size_t>(unshifted_val - min_val);
  if (index >= log_hist.size())
    return log(MIN_PROB);
  return log_hist[index];
}

double 
DiscEmpSplitDistro::log_likelihood(const double unshifted_val,
                                   const double scale) const {
    cerr << "DiscEmpSplitDistro::log_likelihood: "
         << " PLACE HOLDER. NOT WELL TESTED" << endl;

    const size_t index = static_cast<size_t>(unshifted_val / scale - min_val);
  if (index >= log_hist.size())
    return log(MIN_PROB);
  return log_hist[index];
}



DiscEmpSplitDistro::DiscEmpSplitDistro(const DiscEmpSplitDistro &rhs) : 
  SplitDistro_(rhs.params), 
  log_hist(rhs.log_hist),
  hist(rhs.hist),
  cumulative(rhs.cumulative) {
  // TODO: need to set the required parameters based on values in "params" 
}


DiscEmpSplitDistro&
DiscEmpSplitDistro::operator=(const DiscEmpSplitDistro &rhs) {
  if (this != &rhs) {
    this->SplitDistro_::params = rhs.SplitDistro_::params;
    // copy tables
    hist = rhs.hist;
    log_hist = rhs.log_hist;
    cumulative = rhs.cumulative;
  }
  return *this;
}

#include "Smoothing.hpp"

void
DiscEmpSplitDistro::estimate_params_ml(const vector<double> &vals_a,
				       const vector<double> &vals_b) {
  assert(vals_a.size() == vals_b.size());
  
  const size_t lim = vals_a.size();

  vector<double> diffs(vals_a);
  for (size_t i = 0; i < vals_b.size(); ++i)
    diffs[i] -= vals_b[i];

  params[0] = *max_element(diffs.begin(), diffs.end()) + smoothing_parameter;
  params[1] = *min_element(diffs.begin(), diffs.end()) - smoothing_parameter;
  
  n_classes = static_cast<size_t>(params[0] - params[1]) + 1;
  max_val = params[0];
  min_val = params[1];
  params[2] = accumulate(diffs.begin(), diffs.end(), 0.0)/lim;

  vector<double> shifted_vals(diffs);
  for (size_t i = 0; i < shifted_vals.size(); ++i)
    shifted_vals[i] -= min_val;
  
  // build histogram
  make_hist(shifted_vals, n_classes, max_val, min_val, hist);
  
  // make the mids
  vector<double> mids(hist.size());
  for (size_t i = 0; i < mids.size(); ++i)
    mids[i] = i;// + 0.5;
  
  // smooth histogram
  vector<double> smooth_hist;
  LocalLinearRegression(smoothing_parameter, mids, hist, mids, smooth_hist);
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
DiscEmpSplitDistro::estimate_params_ml(const vector<double> &vals_a,
				       const vector<double> &vals_b,
				       const vector<double> &probs) {

  const size_t lim = vals_a.size();
  assert(lim == probs.size() && lim == vals_b.size());

  vector<double> diffs(vals_a);
  for (size_t i = 0; i < vals_b.size(); ++i)
    diffs[i] -= vals_b[i];
  
  params[0] = *max_element(diffs.begin(), diffs.end()) + smoothing_parameter;
  params[1] = *min_element(diffs.begin(), diffs.end()) - smoothing_parameter;

  max_val = params[0];
  min_val = params[1];
  n_classes = static_cast<size_t>(max_val - min_val) + 1;
  params[2] = inner_product(diffs.begin(), diffs.end(), probs.begin(), 0.0)/
    accumulate(probs.begin(), probs.end(), 0.0);
  
  vector<double> shifted_vals(diffs);
  for (size_t i = 0; i < shifted_vals.size(); ++i)
    shifted_vals[i] -= min_val;
  
  // build histogram
  make_weighted_hist(shifted_vals, probs, n_classes, max_val, min_val, hist);
  
  // make the mids
  vector<double> mids(hist.size());
  for (size_t i = 0; i < mids.size(); ++i)
    mids[i] = i;// + 0.5;
  
  // smooth histogram
  vector<double> smooth_hist;
  LocalLinearRegression(smoothing_parameter, mids, hist, mids, smooth_hist);
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
DiscEmpSplitDistro::estimate_params_ml(const vector<double> &vals_a,
                                       const vector<double> &vals_b,
                                       const vector<double> &scales,
                                       const vector<double> &probs) 
{
    cerr << "DiscEmpSplitDistro::estimate_params_ml:  "
         << "Using an untested function" << endl;
    vector<double> v_a(vals_a), v_b(vals_b);

    for (size_t i = 0; i < vals_a.size(); ++i)
    {
        v_a[i] /= scales[i];
        const size_t ca = static_cast<size_t>(v_a[i]);
        v_a[i] = ca + static_cast<size_t>(v_a[i] - ca > 0.5);

        v_b[i] /= scales[i];
        const size_t cb = static_cast<size_t>(v_b[i]);
        v_a[i] = ca + static_cast<size_t>(v_b[i] - cb > 0.5);
    }
    
    estimate_params_ml(v_a, v_b, probs);
}
