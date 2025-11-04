/*
 * Copyright (C) 2012 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

template <class Distro_Type>
void
read_param_file(const std::string &infile, const size_t n,
                std::vector<double> &start_trans,
                std::vector<std::vector<double>> &trans,
                std::vector<double> &end_trans,
                std::vector<Distro_Type> &distros) {

  std::ifstream in(infile);
  if (!in)
    throw std::runtime_error("failed to open: " + infile);

  std::size_t n_from_file{};
  in >> n_from_file;
  if (n_from_file != n)
    throw std::runtime_error("Mismatching number of states");

  std::string tmp_str;
  std::getline(in, tmp_str);

  distros.clear();
  for (size_t i = 0; i < n; ++i) {
    std::string tmp_str;
    std::getline(in, tmp_str);
    distros.push_back(Distro_Type(tmp_str));
  }

  start_trans.resize(n);
  for (size_t i = 0; i < n; ++i)
    in >> start_trans[i];

  trans.resize(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      in >> trans[i][j];

  end_trans.resize(n);
  for (size_t i = 0; i < n; ++i)
    in >> end_trans[i];
}

template <class Distro_Type>
void
write_param_file(const std::string &outfile, const size_t &n,
                 const std::vector<double> &start_trans,
                 const std::vector<std::vector<double>> &trans,
                 const std::vector<double> &end_trans,
                 const std::vector<Distro_Type> &distros) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open file: " + outfile);

  out << "# number of states\n";
  out << n << '\n';

  out << "\n# emmission distributions\n";
  std::copy(std::cbegin(distros), std::cend(distros),
            std::ostream_iterator<Distro_Type>(out, "\n"));

  out << "\n# start to state transition probabilities\n";
  std::copy(std::cbegin(start_trans), std::cend(start_trans),
            std::ostream_iterator<double>(out, "\n"));

  out << "\n# state transition probabilities\n";
  for (size_t i = 0; i < n; ++i) {
    std::copy(std::cbegin(trans[i]), std::cend(trans[i]),
              std::ostream_iterator<double>(out, "\t"));
    out << '\n';
  }

  out << "\n# states to end transition probabilities\n";
  std::copy(std::cbegin(end_trans), std::cend(end_trans),
            std::ostream_iterator<double>(out, "\n"));
}
