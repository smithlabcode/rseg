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

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cmath>
#include <numeric>

#include "SelectBinSize.hpp"
#include "ReadCounts.hpp"
#include "rseg_utils.hpp"

using std::cerr;
using std::endl;
using std::vector;

size_t
select_bin_size(const size_t n_reads, const size_t genome_size,
		const double alpha, const double theta) {
  
  const size_t r = n_reads;
  const size_t G = genome_size;
  
  size_t b_max = 20000;

  size_t b_high = b_max;
  size_t b_low = 0;
  size_t b = b_high/2;

  const double half_alpha = alpha/2;

  while (b_high - b_low > 1) {
    size_t n = G/b;
    double chi_low = gsl_cdf_chisq_Pinv(half_alpha, 2.0/n*r)/(2*r);
    double chi_high = gsl_cdf_chisq_Qinv(half_alpha, 2.0/n*r + 2.0)/(2*r);
    double score = std::max(chi_high/n, chi_low/n);
    if (score > 1.0/theta) {
      b_low = b;
      b = (b + b_high)/2;
    }
    else {
      b_high = b;
      b = (b_low + b)/2;
    }
  }
  return b;
}

size_t
select_bin_size_naive(const vector<SimpleGenomicRegion> &regions,
                      const vector<vector<SimpleGenomicRegion> > &reads,
                      const vector<vector<SimpleGenomicRegion> > &deads)
{
    
  size_t total_genome_size = 0;
  for (size_t i = 0; i < regions.size(); ++i)
    total_genome_size += regions[i].get_width();
  size_t total_dead_size = 0;
  size_t n_reads = 0;
  for (size_t i = 0; i < deads.size(); ++i) 
    {
      for (size_t j = 0; j < deads[i].size(); ++j)
	total_dead_size += deads[i][j].get_width();
      n_reads += reads[i].size();
    }
    
  const size_t genome_size = total_genome_size - total_dead_size;

  const size_t b = static_cast<size_t>(5.0 * genome_size / n_reads);
    
  return std::max(b, static_cast<size_t>(1));
}


size_t
select_bin_size_waterman(const vector<SimpleGenomicRegion> &regions,
                         const vector<vector<SimpleGenomicRegion> > &reads,
                         const vector<vector<SimpleGenomicRegion> > &deads,
                         const bool smooth)
{
    
  const double ref_genome_size = 2287968180;
  const double ref_read_num = 10000000;
  const double ref_bin_size = smooth ? 1000 : 400;
  const double exponent = smooth ? 1.0/6.0 : 1.0/4.0;

  size_t total_genome_size = 0;
  for (size_t i = 0; i < regions.size(); ++i)
    total_genome_size += regions[i].get_width();
  size_t total_dead_size = 0;
  size_t n_reads = 0;
  for (size_t i = 0; i < deads.size(); ++i) 
    {
      for (size_t j = 0; j < deads[i].size(); ++j)
	total_dead_size += deads[i][j].get_width();
      n_reads += reads[i].size();
    }
  const size_t genome_size = total_genome_size - total_dead_size;

  const size_t b =
    static_cast<size_t>(ref_bin_size /
			pow(n_reads/ref_read_num*ref_genome_size/genome_size,
			    exponent));
  return std::max(b, static_cast<size_t>(1));
}

size_t 
select_bin_size_hideaki(const vector<SimpleGenomicRegion> &regions,
                        const vector<vector<SimpleGenomicRegion> > &reads,
                        const vector<vector<SimpleGenomicRegion> > &deads,
                        const bool smooth)
{

  const double ref_genome_size = 2287968180;
  const double ref_read_num = 10000000;
  const double ref_bin_size = smooth ? 1000 : 400;
  const double exponent = smooth ? 1.0/3.0 : 1.0/2.0;

  size_t total_genome_size = 0;
  for (size_t i = 0; i < regions.size(); ++i)
    total_genome_size += regions[i].get_width();
  size_t total_dead_size = 0;
  size_t n_reads = 0;
  for (size_t i = 0; i < deads.size(); ++i) 
    {
      for (size_t j = 0; j < deads[i].size(); ++j)
	total_dead_size += deads[i][j].get_width();
      n_reads += reads[i].size();
    }
  const size_t genome_size = total_genome_size - total_dead_size;

  const size_t b =
    static_cast<size_t>(ref_bin_size /
			pow(n_reads/ref_read_num*ref_genome_size/genome_size,
			    exponent));
  return std::max(b, static_cast<size_t>(1));
}

void
get_mean_var(const vector<double>  &vals, double &mean, double &var)
{
  mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    
  var = 0.0;
  for (size_t i = 0; i < vals.size(); ++i)
    var += pow(vals[i] - mean, 2.0);
  var /= vals.size();
}


size_t 
select_bin_size_hideaki_emp(const vector<SimpleGenomicRegion> &regions,
                            const vector<vector<SimpleGenomicRegion> > &reads,
                            const vector<vector<SimpleGenomicRegion> > &deads,
                            const double max_dead_proportion)
{

  size_t bin_size_low = 10;
  size_t bin_size_high = 20000;
    
  while (bin_size_low < bin_size_high - 3)
    {
      vector<vector<double> > tmp_read_bins;
      vector<vector<SimpleGenomicRegion> > bin_boundaries;
      vector<double> read_bins;
      vector<size_t> reset_points;
        
      size_t b_one_third(bin_size_low + (bin_size_high - bin_size_low) / 3);
      BinReadsCorrectDeadZones(reads, deads, regions, b_one_third,
			       max_dead_proportion,
			       bin_boundaries, tmp_read_bins);
      collapse_read_bins(tmp_read_bins, read_bins, reset_points);

      double mean_one_third, var_one_third, cost_one_third;
      get_mean_var(read_bins, mean_one_third, var_one_third);
      cost_one_third = (2*mean_one_third - var_one_third)
	/ pow(static_cast<double>(b_one_third), 2.0);

      size_t b_two_third(bin_size_high - (bin_size_high - bin_size_low) / 3);
      BinReadsCorrectDeadZones(reads, deads, regions, b_two_third,
			       max_dead_proportion,
			       bin_boundaries, tmp_read_bins);
      collapse_read_bins(tmp_read_bins, read_bins, reset_points);
        
      double mean_two_third, var_two_third, cost_two_third;
      get_mean_var(read_bins, mean_two_third, var_two_third);
      cost_two_third = (2*mean_two_third - var_two_third)
	/ pow(static_cast<double>(b_two_third), 2.0);

      //         cerr << b_one_third << "\t" << cost_one_third << endl
      //              << b_two_third << "\t" << cost_two_third << endl;

      if (cost_one_third > cost_two_third)
	bin_size_low = b_one_third;
      else if (cost_one_third < cost_two_third)
	bin_size_high = b_two_third;
      else if (fabs(cost_one_third - cost_two_third) < 1e-20)
        {
	  bin_size_low = b_one_third;
	  bin_size_high = b_two_third;
        }
    }
    
  return (bin_size_low + bin_size_high) / 2;
}













