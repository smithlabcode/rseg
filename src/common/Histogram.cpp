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

#include "Histogram.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

using std::string;
using std::cerr;
using std::endl;

void
Histogram::set_scale(const double new_mean) {
  scale = new_mean/gsl_histogram_mean(the_hist);
}


size_t
Histogram::round_frac(const double v) const {
  double int_part = std::floor(v);
  double frac_part = v - int_part;
  if (rng.runif(0.0, 1.0) > frac_part) {
//     cerr << v << "\t" << static_cast<size_t>(int_part + 1) << endl;
    return static_cast<size_t>(int_part + 1);
  }
  else {
//     cerr << v << "\t" << static_cast<size_t>(int_part) << endl;
    return static_cast<size_t>(int_part);
  }
}

double
Histogram::sample() const {
  const double s = scale*gsl_histogram_pdf_sample(sampler, rng.runif(0.0, 1.0));
//   return s;
  double int_part = std::floor(s);
  double frac_part = s - int_part;
  if (rng.runif(0.0,1.0) < frac_part)
    return int_part + 1;
  else
    return int_part;
  // return scale*static_cast<double>(gsl_histogram_pdf_sample(sampler, rng.runif(0.0, 1.0)));
}

size_t
Histogram::hist_file_size(string hist_file) {
  std::ifstream in(hist_file.c_str(), std::ios::binary);
  char c = 0;
  size_t hist_size = 0;
  while (in.read(&c, 1))
    if (c == '\n')
      hist_size++;
  in.close();
  return hist_size;
}


Histogram::Histogram(const string hist_file, size_t seed) : 
  scale(1.0), the_hist(0), sampler(0), rng(seed) {
  hist_size = hist_file_size(hist_file);
  the_hist = gsl_histogram_alloc(hist_size);
  FILE *f = fopen(hist_file.c_str(), "r");
  if (gsl_histogram_fscanf(f, the_hist)) {
    cerr << "could not read histogram: " << hist_file << endl;
    exit(EXIT_FAILURE);
  }
  fclose(f);
  
  sampler = gsl_histogram_pdf_alloc(hist_size);
  gsl_histogram_pdf_init(sampler, the_hist);
}
