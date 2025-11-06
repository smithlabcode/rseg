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

#ifndef TWO_STATE_SCALE_HMM_HPP
#define TWO_STATE_SCALE_HMM_HPP

#include <utility>
#include <vector>

struct Distro;

struct TwoStateScaleHMM {
  auto
  ViterbiDecoding(const std::vector<double> &values,
                  const std::vector<double> &scales,
                  const std::vector<std::size_t> &reset_points,
                  const std::vector<double> &start_trans,
                  const std::vector<std::vector<double>> &trans,
                  const std::vector<double> &end_trans, const Distro &fg_distro,
                  const Distro &bg_distro,
                  std::vector<bool> &ml_classes) const -> double;

  auto
  BaumWelchTraining(const std::vector<double> &values,
                    const std::vector<double> &scales,
                    const std::vector<std::size_t> &reset_points,
                    std::vector<double> &start_trans,
                    std::vector<std::vector<double>> &trans,
                    std::vector<double> &end_trans, Distro &fg_distro,
                    Distro &bg_distro) const -> double;

  auto
  PosteriorDecoding(const std::vector<double> &values,
                    const std::vector<double> &scales,
                    const std::vector<std::size_t> &reset_points,
                    const std::vector<double> &start_trans,
                    const std::vector<std::vector<double>> &trans,
                    const std::vector<double> &end_trans,
                    const Distro &fg_distro, const Distro &bg_distro,
                    std::vector<bool> &classes,
                    std::vector<double> &llr_scores) const -> double;

  [[nodiscard]] auto
  PosteriorScores(
    const std::vector<double> &values, const std::vector<double> &scales,
    const std::vector<std::size_t> &reset_points,
    const std::vector<double> &start_trans,
    const std::vector<std::vector<double>> &trans,
    const std::vector<double> &end_trans, const Distro &fg_distro,
    const Distro &bg_distro,
    const std::vector<bool> &classes) const -> std::vector<double>;

  [[nodiscard]] auto
  PosteriorScores(const std::vector<double> &values,
                  const std::vector<double> &scales,
                  const std::vector<std::size_t> &reset_points,
                  const std::vector<double> &start_trans,
                  const std::vector<std::vector<double>> &trans,
                  const std::vector<double> &end_trans, const Distro &fg_distro,
                  const Distro &bg_distro,
                  const bool class_id) const -> std::vector<double>;

  [[nodiscard]] auto
  TransitionPosteriors(
    const std::vector<double> &values, const std::vector<double> &scales,
    const std::vector<std::size_t> &reset_points,
    const std::vector<double> &start_trans,
    const std::vector<std::vector<double>> &trans,
    const std::vector<double> &end_trans, const Distro &fg_distro,
    const Distro &bg_distro,
    const std::size_t transition) const -> std::vector<double>;

  static constexpr std::size_t FG_TO_FG_TRANSITION = 0;
  static constexpr std::size_t FG_TO_BG_TRANSITION = 1;
  static constexpr std::size_t BG_TO_FG_TRANSITION = 2;
  static constexpr std::size_t BG_TO_BG_TRANSITION = 3;

  // clang-format off
  double
  ViterbiDecoding(const std::vector<double> &values,
                  const std::vector<double> &scales,
                  const std::vector<size_t> &reset_points,
                  const double p_sf, const double p_sb,
                  const double p_ff, const double p_fb,
                  const double p_ft, const double p_bf,
                  const double p_bb, const double p_bt,
                  const Distro &fg_distro,
                  const Distro &bg_distro,
                  std::vector<bool> &ml_classes) const;
  // clang-format on

  double
  BaumWelchTraining(const std::vector<double> &values,
                    const std::vector<double> &scales,
                    const std::vector<std::size_t> &reset_points, double &p_sf,
                    double &p_sb, double &p_ff, double &p_fb, double &p_ft,
                    double &p_bf, double &p_bb, double &p_bt, Distro &fg_distro,
                    Distro &bg_distro) const;

  double
  PosteriorDecoding(const std::vector<double> &values,
                    const std::vector<double> &scales,
                    const std::vector<std::size_t> &reset_points, double p_sf,
                    double p_sb, double p_ff, double p_fb, double p_ft,
                    double p_bf, double p_bb, double p_bt,
                    const Distro &fg_distro, const Distro &bg_distro,
                    std::vector<bool> &classes,
                    std::vector<double> &llr_scores) const;

  [[nodiscard]] auto
  PosteriorScores(
    const std::vector<double> &values, const std::vector<double> &scales,
    const std::vector<std::size_t> &reset_points, double p_sf, double p_sb,
    double p_ff, double p_fb, double p_ft, double p_bf, double p_bb,
    double p_bt, const Distro &fg_distro, const Distro &bg_distro,
    const std::vector<bool> &classes) const -> std::vector<double>;

  [[nodiscard]] auto
  PosteriorScores(const std::vector<double> &values,
                  const std::vector<double> &scales,
                  const std::vector<std::size_t> &reset_points, double p_sf,
                  double p_sb, double p_ff, double p_fb, double p_ft,
                  double p_bf, double p_bb, double p_bt,
                  const Distro &fg_distro, const Distro &bg_distro,
                  const bool class_id) const -> std::vector<double>;

  [[nodiscard]] auto
  TransitionPosteriors(
    const std::vector<double> &values, const std::vector<double> &scales,
    const std::vector<std::size_t> &reset_points, double p_sf, double p_sb,
    double p_ff, double p_fb, double p_ft, double p_bf, double p_bb,
    double p_bt, const Distro &fg_distro, const Distro &bg_distro,
    const std::size_t transition) const -> std::vector<double>;

  [[nodiscard]] auto
  single_iteration(const std::vector<double> &values,
                   const std::vector<double> &scales,
                   const std::vector<std::size_t> &reset_points,
                   std::vector<std::pair<double, double>> &forward,
                   std::vector<std::pair<double, double>> &backward,
                   double &p_sf, double &p_sb, double &p_ff, double &p_fb,
                   double &p_ft, double &p_bf, double &p_bb, double &p_bt,
                   Distro &fg_distro, Distro &bg_distro) const -> double;

  [[nodiscard]] auto
  forward_algorithm(const std::vector<double> &vals,
                    const std::vector<double> &scales, const std::size_t start,
                    const std::size_t end, const double lp_sf,
                    const double lp_sb, const double lp_ff, const double lp_fb,
                    const double lp_ft, const double lp_bf, const double lp_bb,
                    const double lp_bt, const Distro &fg_distro,
                    const Distro &bg_distro,
                    std::vector<std::pair<double, double>> &f) const -> double;

  [[nodiscard]] auto
  backward_algorithm(const std::vector<double> &vals,
                     const std::vector<double> &scales, const std::size_t start,
                     const std::size_t end, const double lp_sf,
                     const double lp_sb, const double lp_ff, const double lp_fb,
                     const double lp_ft, const double lp_bf, const double lp_bb,
                     const double lp_bt, const Distro &fg_distro,
                     const Distro &bg_distro,
                     std::vector<std::pair<double, double>> &b) const -> double;

  auto
  estimate_emissions(const std::vector<std::pair<double, double>> &f,
                     const std::vector<std::pair<double, double>> &b,
                     std::vector<double> &fg_probs,
                     std::vector<double> &bg_probs) const -> void;

  auto
  estimate_transitions(
    const std::vector<double> &vals, const std::vector<double> &scales,
    const std::size_t start, const std::size_t end,
    const std::vector<std::pair<double, double>> &f,
    const std::vector<std::pair<double, double>> &b, const double total,
    const Distro &fg_distro, const Distro &bg_distro, const double lp_ff,
    const double lp_fb, const double lp_bf, const double lp_bb,
    const double lp_ft, const double lp_bt, std::vector<double> &ff_vals,
    std::vector<double> &fb_vals, std::vector<double> &bf_vals,
    std::vector<double> &bb_vals) const -> void;

  double MIN_PROB{};
  double tolerance{};
  std::size_t max_iterations{};
  bool VERBOSE{};
};

#endif
