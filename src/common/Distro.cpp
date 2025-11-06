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
#include "lgamma.hpp"
#include "smithlab_utils.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

// DISTRO__

double
Distro_::log_sum_log_vec(const std::vector<double> &vals, std::size_t limit) {
  const auto x = std::max_element(std::cbegin(vals), std::cend(vals) + limit);
  const double max_val = *x;
  const std::size_t max_idx = x - std::cbegin(vals);
  double sum = 1.0;
  for (auto i = 0u; i < limit; ++i) {
    if (i != max_idx) {
      sum += std::exp(vals[i] - max_val);
      assert(std::isfinite(sum));
    }
  }
  return max_val + std::log(sum);
}

Distro_::Distro_() {}

Distro_::Distro_(const std::vector<double> p) : params(p) {}

Distro_::~Distro_() {}

std::string
Distro_::tostring() const {
  std::ostringstream os;
  if (!params.empty())
    os << std::setprecision(4) << params.front();
  for (std::size_t i = 1; i < std::size(params); ++i)
    os << " " << std::setprecision(4) << params[i];
  return os.str();
}

double
Distro_::log_likelihood(std::vector<double>::const_iterator a,
                        std::vector<double>::const_iterator b) const {
  double l = 0;
  for (; a < b; ++a)
    l += log_likelihood(*a);
  return l;
}

double
Distro_::operator()(const double val) const {
  return std::exp(log_likelihood(val));
}

double
Distro_::operator()(const std::vector<double> &vals) const {
  const std::size_t lim = std::size(vals);
  double l = 1;
  for (std::size_t i = 0; i < lim; ++i)
    l *= operator()(vals[i]);
  return l;
}

double
Distro_::log_likelihood(const std::vector<double> &vals) const {
  double l = 0;
  const std::size_t lim = std::size(vals);
  for (std::size_t i = 0; i < lim; ++i)
    l += log_likelihood(vals[i]);
  return l;
}

double
Distro_::log_likelihood(const std::vector<double> &vals,
                        const std::vector<double> &scales) const {
  double l = 0;
  const std::size_t lim = std::size(vals);
  for (std::size_t i = 0; i < lim; ++i)
    l += log_likelihood(vals[i], scales[i]);
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
Distro::has_params(const std::string &name) {
  return (std::size(smithlab::split(name, ",")) >= 2);
}

Distro::Distro(const std::string &n, const std::string &params) :
  name(n), d(distro_factory(n, params)) {}

Distro::Distro(const std::string &n, const std::vector<double> &params) :
  name(n), d(distro_factory(name)) {
  d->set_params(params);
}

Distro::Distro(const std::string &s) : d(distro_factory(s)) {
  std::vector<std::string> name_split;
  if (s.find(",") == std::string::npos)  // whitespaces seperated
    name_split = smithlab::split_whitespace_quoted(s);
  else  // comma seperated
    name_split = smithlab::split(s, ",");
  name = name_split[0];
  name = name_split[0];
}

Distro::Distro(const Distro &rhs) :
  name(rhs.name), d(distro_factory(rhs.name)) {
  std::vector<double> tmp_params(rhs.get_params());
  d->set_params(tmp_params);
}

Distro &
Distro::operator=(const Distro &rhs) {
  if (this != &rhs) {
    name = rhs.name;
    d = distro_factory(rhs.name);
    std::vector<double> tmp_params(rhs.get_params());
    d->set_params(tmp_params);
  }
  return *this;
}

Distro::~Distro() {
  if (d)
    delete d;
}

double
Distro::operator()(double val) const {
  return (*d)(val);
}

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
Distro::log_likelihood(double val) const {
  return d->log_likelihood(val);
}

double
Distro::log_likelihood(const double &val, const double &scale) const {
  return d->log_likelihood(val, scale);
}

std::string
Distro::tostring() const {
  return name + std::string(" ") + d->tostring();
}

std::ostream &
operator<<(std::ostream &s, const Distro &distro) {
  return s << distro.tostring();
}

double
Distro::log_sum_log_vec(const std::vector<double> &vals, std::size_t limit) {
  return Distro_::log_sum_log_vec(vals, limit);
}

////////////////////////////////////////////////////////////////////////
//
// DISTRO FACTORY

Distro_ *
distro_factory(std::string name, std::string params) {
  Distro_ *distro;
  if (name == "std::exp")
    distro = new Std::ExpDistro();
  else if (name == "pois")
    distro = new PoisDistro();
  else if (name == "nbd")
    distro = new NegBinomDistro();
  else
    throw std::runtime_error("bad distribution name: " + name);

  std::vector<std::string> params_split = smithlab::split(params, ",");
  if (std::size(params_split) != distro->required_params())
    throw std::runtime_error(
      "bad number of params: " + std::to_string(std::size(params_split)) +
      " for distro: " + name);
  else {
    std::vector<double> params_vec;
    for (std::size_t i = 0; i < std::size(params_split); ++i)
      params_vec.push_back(atof(params_split[i].c_str()));
    distro->set_params(params_vec);
  }
  return distro;
}

Distro_ *
distro_factory(std::string name_arg) {
  std::vector<std::string> name_split;
  if (name_arg.find(",") == std::string::npos)  // whitespaces seperated
    name_split = smithlab::split_whitespace_quoted(name_arg);
  else  // comma seperated
    name_split = smithlab::split(name_arg, ",");

  const std::string name = name_split.front();

  Distro_ *distro;
  if (name == "std::exp")
    distro = new Std::ExpDistro();
  else if (name == "pois")
    distro = new PoisDistro();
  else if (name == "nbd")
    distro = new NegBinomDistro();
  else
    throw std::runtime_error("bad distribution name \"" + name + "\"");

  if (std::size(name_split) > 1) {

    std::vector<std::string> params_split(
      std::vector<std::string>(name_split.begin() + 1, name_split.end()));
    if (std::size(params_split) != distro->required_params())
      throw std::runtime_error(
        "bad number of params: " + std::to_string(std::size(params_split)) +
        " for distro: " + name);
    else {
      std::vector<double> params_vec;
      for (std::size_t i = 0; i < std::size(params_split); ++i)
        params_vec.push_back(atof(params_split[i].c_str()));
      distro->set_params(params_vec);
    }
  }
  return distro;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

double
Std::ExpDistro::log_likelihood(const double val) const {
  return -std::log(params[0]) - val / params[0];
}

double
Std::ExpDistro::log_likelihood(const double &val, const double &scale) const {
  //// TEST NEEDED
  return -std::log(params[0] * scale) - val / params[0] / scale;
}

Std::ExpDistro::Std::ExpDistro(const Std::ExpDistro &rhs) :
  Distro_(rhs.params) {}

Std::ExpDistro &
Std::ExpDistro::operator=(const Std::ExpDistro &rhs) {
  if (this != &rhs) {
    Distro_::params = rhs.Distro_::params;
  }
  return *this;
}

void
Std::ExpDistro::estimate_params_ml(const std::vector<double> &vals) {
  params.front() =
    std::accumulate(vals.begin(), vals.end(), 0.0) / std::size(vals);
}

void
Std::ExpDistro::estimate_params_ml(const std::vector<double> &vals,
                                   const std::vector<double> &scales,
                                   const std::vector<double> &probs) {
  std::vector<double> values(vals.begin(), vals.end());
  for (std::size_t i = 0; i < std::size(values); ++i)
    values[i] /= scales[i];
  if (std::size(probs) == 0)
    estimate_params_ml(vals);
  else
    estimate_params_ml(vals, probs);
}

void
Std::ExpDistro::estimate_params_ml(const std::vector<double> &vals,
                                   const std::vector<double> &probs) {
  const std::size_t lim = std::size(vals);
  if (std::size(workspace_vals) < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (std::size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = std::log(probs[i]);
    workspace_vals[i] = std::log(vals[i]) + std::log(probs[i]);
  }
  const double prob_sum = std::exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = std::exp(log_sum_log_vec(workspace_vals, lim)) / prob_sum;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

double
PoisDistro::log_likelihood(const double val) const {
  return -params.front() + val * std::log(params.front()) - lfact(val);
}

double
PoisDistro::log_likelihood(const double &val, const double &scale) const {
  const double lambda = params[0] * scale;
  return -lambda + val * std::log(lambda) - lfact(val);
}

PoisDistro::PoisDistro(const PoisDistro &rhs) : Distro_(rhs.params) {}

PoisDistro &
PoisDistro::operator=(const PoisDistro &rhs) {
  if (this != &rhs) {
    Distro_::params = rhs.Distro_::params;
  }
  return *this;
}

void
PoisDistro::estimate_params_ml(const std::vector<double> &vals) {
  params.front() =
    std::accumulate(vals.begin(), vals.end(), 0.0) / std::size(vals);
}

void
PoisDistro::estimate_params_ml(const std::vector<double> &vals,
                               const std::vector<double> &probs) {
  const std::size_t lim = std::size(vals);
  if (std::size(workspace_vals) < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (std::size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = std::log(probs[i]);
    workspace_vals[i] = std::log(vals[i]) + std::log(probs[i]);
  }
  const double prob_sum = std::exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = std::exp(log_sum_log_vec(workspace_vals, lim)) / prob_sum;
}

void
PoisDistro::estimate_params_ml(const std::vector<double> &vals,
                               const std::vector<double> &scales,
                               const std::vector<double> &probs) {
  const std::size_t lim = std::size(vals);
  if (std::size(workspace_vals) < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (std::size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = std::log(probs[i]) + std::log(scales[i]);
    workspace_vals[i] = std::log(vals[i]) + std::log(probs[i]);
  }
  const double prob_sum = std::exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = std::exp(log_sum_log_vec(workspace_vals, lim)) / prob_sum;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

const double NegBinomDistro::max_allowed_alpha = 100;

const double NegBinomDistro::min_allowed_alpha = 1e-20;

const double NegBinomDistro::alpha_allowed_error = 1e-10;

void
NegBinomDistro::set_helpers() {
  n_helper = 1 / params[1];
  p_helper = n_helper / (n_helper + params[0]);
  n_log_p_minus_lngamma_n_helper =
    n_helper * std::log(p_helper) - lgamma(n_helper);
  log_q_helper = std::log(1 - p_helper);
  // TODO: should check that these are valid!!!
}

void
NegBinomDistro::set_params(const std::vector<double> &p) {
  Distro_::set_params(p);
  set_helpers();
}

double
NegBinomDistro::log_likelihood(const double val) const {
  const double P = (lgamma(val + n_helper) - lfact(val)) +
                   n_log_p_minus_lngamma_n_helper + val * log_q_helper;
  if (!std::isfinite(P))
    return -40;
  return P;
}

double
NegBinomDistro::log_likelihood(const double &val, const double &scale) const {
  //// TEST NEEDED
  // alpha is not scaled
  const double scaled_p_helper = n_helper / (n_helper + params[0] * scale);
  const double scaled_n_log_p_minus_lngamma_n_helper =
    n_helper * std::log(scaled_p_helper) - lgamma(n_helper);
  const double scaled_log_q_helper = std::log(1 - scaled_p_helper);

  const double P = (lgamma(val + n_helper) - lfact(val)) +
                   scaled_n_log_p_minus_lngamma_n_helper +
                   val * scaled_log_q_helper;
  if (!std::isfinite(P))
    return -40;
  return P;
}

NegBinomDistro::NegBinomDistro(const NegBinomDistro &rhs) :
  Distro_(rhs.params) {
  set_helpers();
}

NegBinomDistro &
NegBinomDistro::operator=(const NegBinomDistro &rhs) {
  if (this != &rhs) {
    Distro_::params = rhs.Distro_::params;
    set_helpers();
  }
  return *this;
}

static inline double
score_fun_first_term(const std::vector<double> &vals_hist, const double alpha) {
  double sum = 0;
  for (std::size_t i = 0; i < std::size(vals_hist); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (std::size_t j = 0; j < i; ++j)
        inner_sum += j / (1 + alpha * j);
      sum += vals_hist[i] * inner_sum;
    }
  return sum;
}

static inline double
alpha_score_function(const std::vector<double> &vals_hist, const double mu,
                     const double alpha, const double vals_count) {
  const double one_plus_alpha_mu = 1 + alpha * mu;
  return (score_fun_first_term(vals_hist, alpha) / vals_count +
          (std::log(one_plus_alpha_mu) / alpha - mu) / alpha);
}

void
NegBinomDistro::estimate_params_ml(const std::vector<double> &vals) {
  // This is the mu
  params.front() =
    std::accumulate(vals.begin(), vals.end(), 0.0) / std::size(vals);

  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.end());
  std::vector<double> vals_hist(static_cast<std::size_t>(max_value) + 1, 0.0);
  for (std::size_t i = 0; i < std::size(vals); ++i)
    ++vals_hist[static_cast<std::size_t>(vals[i])];

  const double vals_count = std::size(vals);

  const double mu = params.front();
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;

  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error &&
         std::fabs((a_high - a_low) / std::max(a_high, a_low)) >
           alpha_allowed_error) {
    a_mid = (a_low + a_high) / 2;
    const double mid_val =
      alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val) / prev_val);
    prev_val = mid_val;
  }
  params[1] = a_mid;

  set_helpers();
  //   const std::size_t lim = vals.size();

  // This is the mu

  // clang-format off

  // params.front() = std::accumulate(vals.begin(), vals.begin() + lim, 0.0) / lim;
  // const double mu = params.front();
  // const double var = gsl_stats_variance_m(&vals.front(), 1, vals.size(), mu);
  // const double r = (mu * mu) / (var - mu);
  // // const double p = r/(r + params[0]);
  // params[1] = max(0.01, 1 / r);
  // set_helpers();

  // clang-format on
}

void
NegBinomDistro::estimate_params_ml(const std::vector<double> &vals,
                                   const std::vector<double> &probs) {
  //   const std::size_t lim = vals.size();
  //   if (workspace_vals.size() < lim) {
  //     workspace_vals.resize(lim);
  //     workspace_probs.resize(lim);
  //   }
  //   for (std::size_t i = 0; i < lim; ++i) {
  //     workspace_probs[i] = std::log(probs[i]);
  //     workspace_vals[i] = std::log(vals[i]) + std::log(probs[i]);
  //   }
  //   const double vals_count = std::exp(log_sum_log_vec(workspace_probs,
  //   lim)); const double mu = std::exp(log_sum_log_vec(workspace_vals,
  //   lim))/vals_count; const double var = gsl_stats_wvariance_m(&vals.front(),
  //   1, &probs.front(), 1,
  //                                       vals.size(), mu);
  //   const double r = (mu*mu)/(var - mu);
  //   // const double p = r/(r + params[0]);
  //   params[0] = mu;
  //   params[1] = max(0.01, 1/r);
  //   set_helpers();
  const std::size_t lim = std::size(vals);
  if (std::size(workspace_vals) < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (std::size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = std::log(probs[i]);  // - centering_value;
    workspace_vals[i] =
      std::log(vals[i]) + std::log(probs[i]);  // - centering_value;
  }

  const double vals_count = std::exp(log_sum_log_vec(workspace_probs, lim));
  params.front() = std::exp(log_sum_log_vec(workspace_vals, lim)) / vals_count;

  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.begin() + lim);
  std::vector<double> vals_hist(static_cast<std::size_t>(max_value) + 1, 0.0);
  for (std::size_t i = 0; i < lim; ++i)
    vals_hist[static_cast<std::size_t>(vals[i])] += probs[i];

  const double mu = params.front();
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;

  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error &&
         std::fabs((a_high - a_low) / std::max(a_high, a_low)) >
           alpha_allowed_error) {
    a_mid = (a_low + a_high) / 2;
    const double mid_val =
      alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;
    diff = std::fabs((prev_val - mid_val) / std::max(mid_val, prev_val));
    prev_val = mid_val;
  }
  params[1] = a_mid;

  set_helpers();
}

static double
llh_derivative_rt_alpha(const std::vector<double> &vals,
                        const std::vector<double> &scales,
                        const std::vector<double> &probs,
                        const std::vector<double> &vals_hist, const double mu,
                        const double alpha) {
  const double first_term = score_fun_first_term(vals_hist, alpha);

  const double mu_times_alpha = mu * alpha;
  const double alpha_inverse = 1 / alpha;
  const double alpha_square_inverse = pow(alpha_inverse, 2.0);

  double second_term = 0;
  for (std::size_t i = 0; i < std::size(vals); ++i) {
    const double one_plus_extra = 1 + scales[i] * mu_times_alpha;
    second_term +=
      probs[i] * (alpha_square_inverse * std::log(one_plus_extra) -
                  scales[i] * mu * (alpha_inverse + vals[i]) / one_plus_extra);
  }

  return first_term + second_term;
}

void
NegBinomDistro::estimate_params_ml(const std::vector<double> &vals,
                                   const std::vector<double> &scales,
                                   const std::vector<double> &probs) {
  const std::size_t lim = std::size(vals);
  if (std::size(workspace_vals) < lim) {
    workspace_vals.resize(lim);
    workspace_probs.resize(lim);
  }
  for (std::size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = std::log(probs[i]) + std::log(scales[i]);
    workspace_vals[i] = std::log(vals[i]) + std::log(probs[i]);
  }

  // this is mu
  const double mu = std::exp(log_sum_log_vec(workspace_vals, lim)) /
                    std::exp(log_sum_log_vec(workspace_probs, lim));

  // Now for the alpha
  const double max_value = *std::max_element(vals.begin(), vals.end());
  std::vector<double> vals_hist(static_cast<std::size_t>(max_value) + 1, 0.0);
  for (std::size_t i = 0; i < std::size(vals); ++i)
    vals_hist[static_cast<std::size_t>(vals[i])] += probs[i];

  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;

  double a_mid = max_allowed_alpha;
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  while (diff > alpha_allowed_error &&
         std::fabs((a_high - a_low) / std::max(a_high, a_low)) >
           alpha_allowed_error) {
    a_mid = (a_low + a_high) / 2;
    const double mid_val =
      llh_derivative_rt_alpha(vals, scales, probs, vals_hist, mu, a_mid);

    if (mid_val < 0)
      a_high = a_mid;
    else
      a_low = a_mid;

    diff = std::fabs((prev_val - mid_val) / prev_val);
    prev_val = mid_val;
  }

  params[0] = mu;
  params[1] = a_mid;

  set_helpers();
}
