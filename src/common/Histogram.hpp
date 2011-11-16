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

#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "RNG.hpp"
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

class Histogram {
public:
  Histogram(const std::string hist_file, 
	    size_t seed = std::numeric_limits<size_t>::max());
  ~Histogram() {
    if (the_hist) {
      gsl_histogram_free(the_hist);
      the_hist = 0;
    }
    if (sampler) {
      gsl_histogram_pdf_free(sampler);
      sampler = 0;
    }
  }
  double sample() const;
  double get_mean() const {return gsl_histogram_mean(the_hist);}
  void set_scale(const double new_mean);
  size_t round_frac(const double v) const;
private:
  size_t hist_size;
  double scale;
  gsl_histogram *the_hist;
  gsl_histogram_pdf *sampler;
  Runif rng;
  static size_t hist_file_size(std::string hist_file);
};

#endif
