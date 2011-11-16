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

#ifndef EVALUATE_BOUNDARIES_HPP
#define EVALUATE_BOUNDARIES_HPP

#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

using std::vector;


struct Domain {
  vector<double> vals;
  size_t state;
  Domain(const vector<double> &v, 
	 size_t start, size_t end, bool s) : state(s) {
    vals.insert(vals.end(), v.begin() + start, v.begin() + end);
  }
  std::string tostring() const;
};


class BoundEval {
public:
  BoundEval(const size_t b, const double bw);
  void 
  evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
	   const vector<vector<bool> > &classes,
	   const vector<vector<double> > &scores,
	   vector<vector<GenomicRegion> > &boundary_scores,
	   vector<vector<GenomicRegion> > &boundary_peaks,
	   vector<vector<size_t> > &boundary_sizes) const;
  void 
  evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
	   const vector<vector<size_t> > &classes,
	   const vector<vector<double> > &scores,
	   vector<vector<GenomicRegion> > &boundary_scores,
	   vector<vector<GenomicRegion> > &boundary_peaks,
	   vector<vector<size_t> > &boundary_sizes) const;
  void 
  evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
	   const vector<vector<bool> > &classes,
	   const vector<vector<double> > &scores,
	   vector<vector<GenomicRegion> > &boundaries) const;

  void
  evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
	   const vector<vector<size_t> > &classes,
	   const vector<vector<double> > &scores,
	   vector<vector<GenomicRegion> > &boundaries) const;

  void  
  evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
           const vector<vector<bool> > &classes,
           const vector<vector<double> > &scores,
           const vector<vector<double> > &fg_to_fg_trans_score,
           const vector<vector<double> > &fg_to_bg_trans_score,
           const vector<vector<double> > &bg_to_fg_trans_score,
           const vector<vector<double> > &bg_to_bg_trans_score,
           vector<vector<GenomicRegion> > &boundaries,
           vector<vector<GenomicRegion> > &boundary_peaks,
           vector<vector<size_t> > &boundary_sizes) const;
    
    void  
    evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
             const vector<size_t> &reset_points,
             const vector<bool> &classes,
             const vector<double> &trans_scores,
             const vector<double> &fg_to_fg_trans_score,
             const vector<double> &fg_to_bg_trans_score,
             const vector<double> &bg_to_fg_trans_score,
             const vector<double> &bg_to_bg_trans_score,
             const double cutoff,
             vector<GenomicRegion> &boundaries) const;

    void  
    evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
             const vector<size_t> &reset_points,
             const vector<size_t> &classes,
             const vector<double> &trans_scores,
             const vector<vector<vector<double> > > &post_trans,
             const double cutoff,
             vector<GenomicRegion> &boundaries) const;
    void  
    evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
             const vector<size_t> &reset_points,
             const vector<bool> &classes,
             const vector<double>  &trans_scores,
             const vector<double> &fg_to_fg_trans_score,
             const vector<double> &fg_to_bg_trans_score,
             const vector<double> &bg_to_fg_trans_score,
             const vector<double> &bg_to_bg_trans_score,
             const double cutoff,
             const bool Both_Domain_Ends,
             vector<GenomicRegion> &boundaries) const;
private:
  size_t boundary_size;
  double bandwidth;
};

#endif
