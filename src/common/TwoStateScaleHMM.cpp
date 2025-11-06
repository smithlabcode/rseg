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

#include "TwoStateScaleHMM.hpp"
#include "Distro.hpp"
#include "log_sum_log.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

double
TwoStateScaleHMM::forward_algorithm(
  const std::vector<double> &vals, const std::vector<double> &scales,
  const size_t start, const size_t end, const double lp_sf, const double lp_sb,
  const double lp_ff, const double lp_fb, const double lp_ft,
  const double lp_bf, const double lp_bb, const double lp_bt,
  const Distro &fg_distro, const Distro &bg_distro,
  std::vector<std::pair<double, double>> &f) const {

  f[start].first = fg_distro.log_likelihood(vals[start], scales[start]) + lp_sf;
  f[start].second =
    bg_distro.log_likelihood(vals[start], scales[start]) + lp_sb;

  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    f[i].first = (fg_distro.log_likelihood(vals[i], scales[i]) +
                  log_sum_log(f[k].first + lp_ff, f[k].second + lp_bf));
    f[i].second = (bg_distro.log_likelihood(vals[i], scales[i]) +
                   log_sum_log(f[k].first + lp_fb, f[k].second + lp_bb));
  }
  return log_sum_log(f[end - 1].first + lp_ft, f[end - 1].second + lp_bt);
}

double
TwoStateScaleHMM::backward_algorithm(
  const std::vector<double> &vals, const std::vector<double> &scales,
  const size_t start, const size_t end, const double lp_sf, const double lp_sb,
  const double lp_ff, const double lp_fb, const double lp_ft,
  const double lp_bf, const double lp_bb, const double lp_bt,
  const Distro &fg_distro, const Distro &bg_distro,
  std::vector<std::pair<double, double>> &b) const {

  b[end - 1].first = lp_ft;
  b[end - 1].second = lp_bt;

  for (size_t k = end - 1; k > start; --k) {
    size_t i = k - 1;
    const double fg_a =
      fg_distro.log_likelihood(vals[k], scales[k]) + b[k].first;
    const double bg_a =
      bg_distro.log_likelihood(vals[k], scales[k]) + b[k].second;
    b[i].first = log_sum_log(fg_a + lp_ff, bg_a + lp_fb);
    b[i].second = log_sum_log(fg_a + lp_bf, bg_a + lp_bb);
  }
  return log_sum_log(
    b[start].first + fg_distro.log_likelihood(vals[start], scales[start]) +
      lp_sf,
    b[start].second + bg_distro.log_likelihood(vals[start], scales[start]) +
      lp_sb);
}

void
TwoStateScaleHMM::estimate_emissions(
  const std::vector<std::pair<double, double>> &f,
  const std::vector<std::pair<double, double>> &b,
  std::vector<double> &fg_probs, std::vector<double> &bg_probs) const {
  for (size_t i = 0; i < b.size(); ++i) {
    const double fg = (f[i].first + b[i].first);
    const double bg = (f[i].second + b[i].second);
    const double denom = log_sum_log(fg, bg);
    fg_probs[i] = std::exp(fg - denom);
    bg_probs[i] = std::exp(bg - denom);
  }
}

void
TwoStateScaleHMM::estimate_transitions(
  const std::vector<double> &vals, const std::vector<double> &scales,
  const size_t start, const size_t end,
  const std::vector<std::pair<double, double>> &f,
  const std::vector<std::pair<double, double>> &b, const double total,
  const Distro &fg_distro, const Distro &bg_distro, const double lp_ff,
  const double lp_fb, const double lp_bf, const double lp_bb,
  [[maybe_unused]] const double lp_ft, [[maybe_unused]] const double lp_bt,
  std::vector<double> &ff_vals, std::vector<double> &fb_vals,
  std::vector<double> &bf_vals, std::vector<double> &bb_vals) const {

  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;

    const double lp_fg = fg_distro.log_likelihood(vals[i], scales[i]) - total;
    const double lp_bg = bg_distro.log_likelihood(vals[i], scales[i]) - total;

    ff_vals[k] = f[k].first + lp_ff + lp_fg + b[i].first;
    fb_vals[k] = f[k].first + lp_fb + lp_bg + b[i].second;

    bf_vals[k] = f[k].second + lp_bf + lp_fg + b[i].first;
    bb_vals[k] = f[k].second + lp_bb + lp_bg + b[i].second;
  }
}

double
TwoStateScaleHMM::single_iteration(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points,
  std::vector<std::pair<double, double>> &forward,
  std::vector<std::pair<double, double>> &backward, double &p_sf, double &p_sb,
  double &p_ff, double &p_fb, double &p_ft, double &p_bf, double &p_bb,
  double &p_bt, Distro &fg_distro, Distro &bg_distro) const {

  std::vector<double> log_fg_expected;
  std::vector<double> log_bg_expected;

  double total_score = 0;

  const double lp_sf = std::log(p_sf);
  const double lp_sb = std::log(p_sb);
  const double lp_ff = std::log(p_ff);
  const double lp_fb = std::log(p_fb);
  const double lp_ft = std::log(p_ft);
  const double lp_bf = std::log(p_bf);
  const double lp_bb = std::log(p_bb);
  const double lp_bt = std::log(p_bt);

  assert(std::isfinite(lp_sf) && std::isfinite(lp_sb) && std::isfinite(lp_ff) &&
         std::isfinite(lp_fb) && std::isfinite(lp_ft) && std::isfinite(lp_bf) &&
         std::isfinite(lp_bb) && std::isfinite(lp_bt));

  // for estimating transitions
  std::vector<double> ff_vals(values.size(), 0);
  std::vector<double> fb_vals(values.size(), 0);
  std::vector<double> bf_vals(values.size(), 0);
  std::vector<double> bb_vals(values.size(), 0);

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, forward);
    [[maybe_unused]] const double backward_score = backward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, backward);

    estimate_transitions(values, scales, reset_points[i], reset_points[i + 1],
                         forward, backward, score, fg_distro, bg_distro, lp_ff,
                         lp_fb, lp_bf, lp_bb, lp_ft, lp_bt, ff_vals, fb_vals,
                         bf_vals, bb_vals);

    total_score += score;
  }

  // Subtracting 1 from the limit of the summation because the final term has
  // no meaning since there is no transition to be counted from the final
  // observation (they all must go to terminal state) SQ: note
  // ff_vals[reset_points[i]] is always euqal to 0 (std::exp() == 1), i.e. the
  // start of each region does not contribute to transition emission estimate
  const double p_ff_new_estimate =
    std::exp(log_sum_log_vec(ff_vals, values.size())) - reset_points.size() + 1;
  const double p_fb_new_estimate =
    std::exp(log_sum_log_vec(fb_vals, values.size())) - reset_points.size() + 1;
  const double p_bf_new_estimate =
    std::exp(log_sum_log_vec(bf_vals, values.size())) - reset_points.size() + 1;
  const double p_bb_new_estimate =
    std::exp(log_sum_log_vec(bb_vals, values.size())) - reset_points.size() + 1;

  // clang-format off
  // const double p_ff_new_estimate =
  //   std::exp(log_sum_log_vec(ff_vals, values.size() - 1));
  // const double p_fb_new_estimate =
  //   std::exp(log_sum_log_vec(fb_vals, values.size() - 1));
  // const double p_bf_new_estimate =
  //   std::exp(log_sum_log_vec(bf_vals, values.size() - 1));
  // const double p_bb_new_estimate =
  //   std::exp(log_sum_log_vec(bb_vals, values.size() - 1));
  // clang-format on

  double denom = (p_ff_new_estimate + p_fb_new_estimate);
  p_ff = p_ff_new_estimate / denom - p_ft / 2.0;
  p_fb = p_fb_new_estimate / denom - p_ft / 2.0;

  if (p_ff < MIN_PROB)
    p_ff = MIN_PROB;

  if (p_fb < MIN_PROB)
    p_fb = MIN_PROB;

  denom = (p_bf_new_estimate + p_bb_new_estimate);
  p_bf = p_bf_new_estimate / denom - p_bt / 2.0;
  p_bb = p_bb_new_estimate / denom - p_bt / 2.0;

  if (p_bf < MIN_PROB)
    p_bf = MIN_PROB;

  if (p_bb < MIN_PROB)
    p_bb = MIN_PROB;

  // for estimating emissions
  std::vector<double> fg_probs(values.size());
  std::vector<double> bg_probs(values.size());
  estimate_emissions(forward, backward, fg_probs, bg_probs);

  fg_distro.estimate_params_ml(values, scales, fg_probs);
  bg_distro.estimate_params_ml(values, scales, bg_probs);

  return total_score;
}

double
TwoStateScaleHMM::BaumWelchTraining(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points, std::vector<double> &start_trans,
  std::vector<std::vector<double>> &trans, std::vector<double> &end_trans,
  Distro &fg_distro, Distro &bg_distro) const {

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  assert(std::all_of(std::cbegin(trans), std::cend(trans),
                     [](auto const &t) { return std::size(t) >= 2; }));

  return BaumWelchTraining(values, scales, reset_points, start_trans[0],
                           start_trans[1], trans[0][0], trans[0][1],
                           end_trans[0], trans[1][0], trans[1][1], end_trans[1],
                           fg_distro, bg_distro);
}

double
TwoStateScaleHMM::BaumWelchTraining(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points, double &p_sf, double &p_sb,
  double &p_ff, double &p_fb, double &p_ft, double &p_bf, double &p_bb,
  double &p_bt, Distro &fg_distro, Distro &bg_distro) const {

  const auto n_values = std::size(values);
  std::vector<std::pair<double, double>> forward(n_values);
  std::vector<std::pair<double, double>> backward(n_values);

  if (VERBOSE)
    std::cout << std::setw(5) << "ITR" << std::setw(10) << "F size"
              << std::setw(10) << "B size" << std::setw(18) << "F PARAMS"
              << std::setw(18) << "B PARAMS" << std::setw(14) << "DELTA"
              << std::endl;

  double prev_total = -std::numeric_limits<double>::max();

  for (size_t i = 0; i < max_iterations; ++i) {
    double p_sf_est = p_sf;
    double p_sb_est = p_sb;
    double p_ff_est = p_ff;
    double p_fb_est = p_fb;
    double p_bf_est = p_bf;
    double p_bb_est = p_bb;
    double p_ft_est = p_ft;
    double p_bt_est = p_bt;

    const auto total =
      single_iteration(values, scales, reset_points, forward, backward,
                       p_sf_est, p_sb_est, p_ff_est, p_fb_est, p_ft_est,
                       p_bf_est, p_bb_est, p_bt_est, fg_distro, bg_distro);

    if ((prev_total - total) / prev_total < tolerance) {
      if (VERBOSE)
        std::cout << "CONVERGED\n" << std::endl;
      break;
    }

    if (VERBOSE) {
      // clang-format off
      std::cout << std::setw(5) << i + 1
                << std::setw(10) << 1 / p_fb_est
                << std::setw(10) << 1 / p_bf_est
                << std::setw(18) << fg_distro.tostring()
                << std::setw(18) << bg_distro.tostring()
                << std::setw(14) << (prev_total - total) / prev_total
                << std::endl;
      // clang-format on
    }

    p_sf = p_sf_est;
    p_sb = p_sb_est;
    p_ff = p_ff_est;
    p_fb = p_fb_est;
    p_bf = p_bf_est;
    p_bb = p_bb_est;
    p_ft = p_ft_est;
    p_bt = p_bt_est;

    prev_total = total;
  }
  return prev_total;
}

[[nodiscard]] auto
TwoStateScaleHMM::PosteriorScores(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points,
  const std::vector<double> &start_trans,
  const std::vector<std::vector<double>> &trans,
  const std::vector<double> &end_trans, const Distro &fg_distro,
  const Distro &bg_distro,
  const std::vector<bool> &classes) const -> std::vector<double> {
  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  assert(std::all_of(std::cbegin(trans), std::cend(trans),
                     [](auto const &t) { return std::size(t) >= 2; }));

  return PosteriorScores(values, scales, reset_points, start_trans[0],
                         start_trans[1], trans[0][0], trans[0][1], end_trans[0],
                         trans[1][0], trans[1][1], end_trans[1], fg_distro,
                         bg_distro, classes);
}

[[nodiscard]] auto
TwoStateScaleHMM::PosteriorScores(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points, double p_sf, double p_sb,
  double p_ff, double p_fb, double p_ft, double p_bf, double p_bb, double p_bt,
  const Distro &fg_distro, const Distro &bg_distro,
  const std::vector<bool> &classes) const -> std::vector<double> {

  double total_score = 0;

  const double lp_sf = std::log(p_sf);
  const double lp_sb = std::log(p_sb);
  const double lp_ff = std::log(p_ff);
  const double lp_fb = std::log(p_fb);
  const double lp_ft = std::log(p_ft);
  const double lp_bf = std::log(p_bf);
  const double lp_bb = std::log(p_bb);
  const double lp_bt = std::log(p_bt);

  assert(std::isfinite(lp_sf) && std::isfinite(lp_sb) && std::isfinite(lp_ff) &&
         std::isfinite(lp_fb) && std::isfinite(lp_ft) && std::isfinite(lp_bf) &&
         std::isfinite(lp_bb) && std::isfinite(lp_bt));

  std::vector<std::pair<double, double>> forward(
    values.size(), std::pair<double, double>(0, 0));
  std::vector<std::pair<double, double>> backward(
    values.size(), std::pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, forward);

    [[maybe_unused]] const double backward_score = backward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, backward);

    total_score += score;
  }

  std::vector<double> llr_scores(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    const double denom = log_sum_log(fg_state, bg_state);
    llr_scores[i] =
      classes[i] ? std::exp(fg_state - denom) : std::exp(bg_state - denom);
  }
  return llr_scores;
}

[[nodiscard]] auto
TwoStateScaleHMM::PosteriorScores(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points,
  const std::vector<double> &start_trans,
  const std::vector<std::vector<double>> &trans,
  const std::vector<double> &end_trans, const Distro &fg_distro,
  const Distro &bg_distro, const bool fg_class) const -> std::vector<double> {
  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  assert(std::all_of(std::cbegin(trans), std::cend(trans),
                     [](auto const &t) { return std::size(t) >= 2; }));

  return PosteriorScores(values, scales, reset_points, start_trans[0],
                         start_trans[1], trans[0][0], trans[0][1], end_trans[0],
                         trans[1][0], trans[1][1], end_trans[1], fg_distro,
                         bg_distro, fg_class);
}

[[nodiscard]] auto
TwoStateScaleHMM::PosteriorScores(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points, double p_sf, double p_sb,
  double p_ff, double p_fb, double p_ft, double p_bf, double p_bb, double p_bt,
  const Distro &fg_distro, const Distro &bg_distro,
  const bool fg_class) const -> std::vector<double> {

  double total_score = 0;

  const double lp_sf = std::log(p_sf);
  const double lp_sb = std::log(p_sb);
  const double lp_ff = std::log(p_ff);
  const double lp_fb = std::log(p_fb);
  const double lp_ft = std::log(p_ft);
  const double lp_bf = std::log(p_bf);
  const double lp_bb = std::log(p_bb);
  const double lp_bt = std::log(p_bt);

  assert(std::isfinite(lp_sf) && std::isfinite(lp_sb) && std::isfinite(lp_ff) &&
         std::isfinite(lp_fb) && std::isfinite(lp_ft) && std::isfinite(lp_bf) &&
         std::isfinite(lp_bb) && std::isfinite(lp_bt));

  std::vector<std::pair<double, double>> forward(
    values.size(), std::pair<double, double>(0, 0));
  std::vector<std::pair<double, double>> backward(
    values.size(), std::pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, forward);

    [[maybe_unused]] const double backward_score = backward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, backward);

    total_score += score;
  }

  std::vector<double> llr_scores(std::size(values));
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    const double denom = log_sum_log(fg_state, bg_state);
    llr_scores[i] =
      fg_class ? std::exp(fg_state - denom) : std::exp(bg_state - denom);
  }
  return llr_scores;
}

[[nodiscard]] auto
TwoStateScaleHMM::TransitionPosteriors(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points,
  const std::vector<double> &start_trans,
  const std::vector<std::vector<double>> &trans,
  const std::vector<double> &end_trans, const Distro &fg_distro,
  const Distro &bg_distro,
  const size_t transition) const -> std::vector<double> {
  assert(std::size(start_trans) >= 2);
  assert(std::size(end_trans) >= 2);
  assert(std::size(trans) >= 2);
  assert(std::all_of(std::cbegin(trans), std::cend(trans),
                     [](auto const &t) { return std::size(t) >= 2; }));
  return TransitionPosteriors(values, scales, reset_points, start_trans[0],
                              start_trans[1], trans[0][0], trans[0][1],
                              end_trans[0], trans[1][0], trans[1][1],
                              end_trans[1], fg_distro, bg_distro, transition);
}

[[nodiscard]] auto
TwoStateScaleHMM::TransitionPosteriors(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points, double p_sf, double p_sb,
  double p_ff, double p_fb, double p_ft, double p_bf, double p_bb, double p_bt,
  const Distro &fg_distro, const Distro &bg_distro,
  const size_t transition) const -> std::vector<double> {
  double total_score = 0;

  const double lp_sf = std::log(p_sf);
  const double lp_sb = std::log(p_sb);
  const double lp_ff = std::log(p_ff);
  const double lp_fb = std::log(p_fb);
  const double lp_ft = std::log(p_ft);
  const double lp_bf = std::log(p_bf);
  const double lp_bb = std::log(p_bb);
  const double lp_bt = std::log(p_bt);

  assert(std::isfinite(lp_sf) && std::isfinite(lp_sb) && std::isfinite(lp_ff) &&
         std::isfinite(lp_fb) && std::isfinite(lp_ft) && std::isfinite(lp_bf) &&
         std::isfinite(lp_bb) && std::isfinite(lp_bt));

  std::vector<std::pair<double, double>> forward(
    values.size(), std::pair<double, double>(0, 0));
  std::vector<std::pair<double, double>> backward(
    values.size(), std::pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, forward);

    [[maybe_unused]] const double backward_score = backward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, backward);

    total_score += score;
  }

  const auto n_values = std::size(values);

  std::vector<double> scores(n_values);
  size_t j = 0;
  for (size_t i = 0; i < n_values; ++i) {
    if (i == reset_points[j]) {
      ++j;
      scores[i] = 0;
    }
    else {
      const double fg_to_fg_state =
        forward[i - 1].first + lp_ff +  // transition
                                        // emission for value i + 1
        fg_distro.log_likelihood(values[i], scales[i]) + backward[i].first;
      const double fg_to_bg_state =
        forward[i - 1].first + lp_fb +
        bg_distro.log_likelihood(values[i], scales[i]) + backward[i].second;
      const double bg_to_fg_state =
        forward[i - 1].second + lp_bf +
        fg_distro.log_likelihood(values[i], scales[i]) + backward[i].first;
      const double bg_to_bg_state =
        forward[i - 1].second + lp_bb +
        bg_distro.log_likelihood(values[i], scales[i]) + backward[i].second;
      const double denom =
        log_sum_log(log_sum_log(fg_to_fg_state, fg_to_bg_state),
                    log_sum_log(bg_to_fg_state, bg_to_bg_state));
      double numerator = fg_to_fg_state;
      if (transition == 1)
        numerator = fg_to_bg_state;
      if (transition == 2)
        numerator = bg_to_fg_state;
      if (transition == 3)
        numerator = bg_to_bg_state;
      scores[i] = std::exp(numerator - denom);
    }
  }
  return scores;
}

double
TwoStateScaleHMM::PosteriorDecoding(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points,
  const std::vector<double> &start_trans,
  const std::vector<std::vector<double>> &trans,
  const std::vector<double> &end_trans, const Distro &fg_distro,
  const Distro &bg_distro, std::vector<bool> &classes,
  std::vector<double> &llr_scores) const {

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  assert(std::all_of(std::cbegin(trans), std::cend(trans),
                     [](auto const &t) { return std::size(t) >= 2; }));

  return PosteriorDecoding(values, scales, reset_points, start_trans[0],
                           start_trans[1], trans[0][0], trans[0][1],
                           end_trans[0], trans[1][0], trans[1][1], end_trans[1],
                           fg_distro, bg_distro, classes, llr_scores);
}

double
TwoStateScaleHMM::PosteriorDecoding(
  const std::vector<double> &values, const std::vector<double> &scales,
  const std::vector<size_t> &reset_points, double p_sf, double p_sb,
  double p_ff, double p_fb, double p_ft, double p_bf, double p_bb, double p_bt,
  const Distro &fg_distro, const Distro &bg_distro, std::vector<bool> &classes,
  std::vector<double> &llr_scores) const {

  double total_score = 0;

  const double lp_sf = std::log(p_sf);
  const double lp_sb = std::log(p_sb);
  const double lp_ff = std::log(p_ff);
  const double lp_fb = std::log(p_fb);
  const double lp_ft = std::log(p_ft);
  const double lp_bf = std::log(p_bf);
  const double lp_bb = std::log(p_bb);
  const double lp_bt = std::log(p_bt);

  assert(std::isfinite(lp_sf) && std::isfinite(lp_sb) && std::isfinite(lp_ff) &&
         std::isfinite(lp_fb) && std::isfinite(lp_ft) && std::isfinite(lp_bf) &&
         std::isfinite(lp_bb) && std::isfinite(lp_bt));

  const auto n_values = std::size(values);
  std::vector<std::pair<double, double>> forward(n_values);
  std::vector<std::pair<double, double>> backward(n_values);

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, forward);

    [[maybe_unused]] const double backward_score = backward_algorithm(
      values, scales, reset_points[i], reset_points[i + 1], lp_sf, lp_sb, lp_ff,
      lp_fb, lp_ft, lp_bf, lp_bb, lp_bt, fg_distro, bg_distro, backward);

    total_score += score;
  }

  classes.resize(n_values);
  llr_scores.resize(n_values);
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    const double denom = log_sum_log(fg_state, bg_state);
    classes[i] = fg_state > bg_state;
    llr_scores[i] =
      classes[i] ? std::exp(fg_state - denom) : std::exp(bg_state - denom);
  }

  return total_score;
}

// Functions for Viterbi training and decoding.

double
TwoStateScaleHMM::ViterbiDecoding(const std::vector<double> &values,
                                  const std::vector<double> &scales,
                                  const std::vector<size_t> &reset_points,
                                  const std::vector<double> &start_trans,
                                  const std::vector<std::vector<double>> &trans,
                                  const std::vector<double> &end_trans,
                                  const Distro &fg_distro,
                                  const Distro &bg_distro,
                                  std::vector<bool> &classes) const {
  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  assert(std::all_of(std::cbegin(trans), std::cend(trans),
                     [](auto const &t) { return std::size(t) >= 2; }));

  return ViterbiDecoding(values, scales, reset_points, start_trans[0],
                         start_trans[1], trans[0][0], trans[0][1], end_trans[0],
                         trans[1][0], trans[1][1], end_trans[1], fg_distro,
                         bg_distro, classes);
}

double
TwoStateScaleHMM::ViterbiDecoding(
  // clang-format off
  const std::vector<double> &values,
  const std::vector<double> &scales,
  const std::vector<size_t> &reset_points,
  const double p_sf, const double p_sb,
  const double p_ff, const double p_fb,
  const double p_ft, const double p_bf,
  const double p_bb, const double p_bt,
  const Distro &fg_distro,
  const Distro &bg_distro,
  std::vector<bool> &ml_classes
  // clang-format on
) const {

  const double lp_sf = std::log(p_sf);
  const double lp_sb = std::log(p_sb);
  const double lp_ff = std::log(p_ff);
  const double lp_fb = std::log(p_fb);
  const double lp_ft = std::log(p_ft);
  const double lp_bf = std::log(p_bf);
  const double lp_bb = std::log(p_bb);
  const double lp_bt = std::log(p_bt);

  // ml_classes = std::vector<bool>(values.size());
  double total = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

    const size_t start = reset_points[i];
    const size_t lim = reset_points[i + 1] - start;

    std::vector<std::pair<double, double>> v(lim,
                                             std::pair<double, double>(0, 0));
    std::vector<std::pair<size_t, size_t>> trace(
      lim, std::pair<size_t, size_t>(0, 0));

    v.front().first =
      lp_sf + fg_distro.log_likelihood(values[start], scales[start]);
    v.front().second =
      lp_sb + bg_distro.log_likelihood(values[start], scales[start]);

    for (size_t j = 1; j < lim; ++j) {

      const double ff = v[j - 1].first + lp_ff;
      const double bf = v[j - 1].second + lp_bf;
      const double fg_log_emmit =
        fg_distro.log_likelihood(values[start + j], scales[start + j]);
      if (ff > bf) {
        v[j].first = fg_log_emmit + ff;
        trace[j].first = 0;
      }
      else {
        v[j].first = fg_log_emmit + bf;
        trace[j].first = 1;
      }

      const double fb = v[j - 1].first + lp_fb;
      const double bb = v[j - 1].second + lp_bb;
      const double bg_log_emmit =
        bg_distro.log_likelihood(values[start + j], scales[start + j]);
      if (fb > bb) {
        v[j].second = bg_log_emmit + fb;
        trace[j].second = 0;
      }
      else {
        v[j].second = bg_log_emmit + bb;
        trace[j].second = 1;
      }
    }
    v.back().first += lp_ft;
    v.back().second += lp_bt;

    std::vector<bool> inner_ml_classes;

    // do the traceback
    size_t prev = 0;
    if (v.back().first > v.back().second) {
      inner_ml_classes.push_back(true);
      prev = trace.back().first;
    }
    else {
      inner_ml_classes.push_back(false);
      prev = trace.back().second;
    }

    for (size_t j = trace.size() - 1; j > 0; --j) {
      const size_t k = j - 1;
      if (prev == 0) {
        inner_ml_classes.push_back(true);
        prev = trace[k].first;
      }
      else {
        inner_ml_classes.push_back(false);
        prev = trace[k].second;
      }
    }

    reverse(inner_ml_classes.begin(), inner_ml_classes.end());
    ml_classes.insert(ml_classes.end(), inner_ml_classes.begin(),
                      inner_ml_classes.end());

    total += std::max(v.back().first, v.back().second);
  }

  return total;
}
