/* MIT License
 *
 * Copyright (c) 2025 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef INTERVAL_HPP_
#define INTERVAL_HPP_

#include <cstdint>
#include <format>
#include <iterator>  // std::size
#include <stdexcept>
#include <string>
#include <vector>

struct Interval {
  std::string chrom;
  std::uint32_t start{};
  std::uint32_t stop{};

  Interval() = default;
  Interval(const std::string &chrom, const std::uint32_t start,
           const std::uint32_t stop) :
    chrom{chrom},
    start{start}, stop{stop} {}

  explicit Interval(const std::string &line) {
    if (!initialize(line.data(), line.data() + std::size(line)))
      throw std::runtime_error("bad interval line: " + line);
  }
  auto
  initialize(const char *, const char *) -> bool;
  auto
  operator<=>(const Interval &) const = default;
};

template <> struct std::formatter<Interval> : std::formatter<std::string> {
  auto
  format(const Interval &i, format_context &ctx) const {
    static constexpr auto fmt = "{}\t{}\t{}";
    return std::formatter<std::string>::format(
      std::format(fmt, i.chrom, i.start, i.stop), ctx);
  }
};

[[nodiscard]] inline auto
size(const Interval &x) {
  return x.stop > x.start ? x.stop - x.start : 0ul;
}

[[nodiscard]] auto
read_intervals(const std::string &intervals_file) -> std::vector<Interval>;

#endif  // INTERVAL_HPP_
