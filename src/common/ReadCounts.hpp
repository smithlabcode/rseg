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

#ifndef READ_COUNTS_HPP
#define READ_COUNTS_HPP

#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

using std::vector;

void
BinReads(const std::vector<std::vector<SimpleGenomicRegion> > &reads,
	 const std::vector<SimpleGenomicRegion> &regions,
	 const size_t bin_size,
	 std::vector<std::vector<double> > &bins);


void
BinReadsExtendOverDeadZones(const std::vector<std::vector<SimpleGenomicRegion> > 
			    &reads, 
			    const std::vector<std::vector<SimpleGenomicRegion> > 
			    &dead_zones,
			    const std::vector<SimpleGenomicRegion> &regions,
			    const size_t bin_size,
			    std::vector<std::vector<SimpleGenomicRegion> > 
			    &boundaries,
			    std::vector<std::vector<double> > &bins);



void
BinReadsCorrectDeadZones(const std::vector<std::vector<SimpleGenomicRegion> > 
			 &reads, 
			 const std::vector<std::vector<SimpleGenomicRegion> > 
			 &dead_zones,
			 const std::vector<SimpleGenomicRegion> &regions,
			 const size_t bin_size,
			 const double max_dead_proportion,
			 std::vector<std::vector<SimpleGenomicRegion> > 
			 &boundaries,
			 std::vector<std::vector<double> > &bins);

void
BinReadsCorrectDeadZones(
    const vector<vector<SimpleGenomicRegion> > &reads, 
    const vector<vector<SimpleGenomicRegion> > &dead_zones,
    const vector<SimpleGenomicRegion> &regions,
    const size_t bin_size,
    const double max_dead_proportion,
    vector<vector<SimpleGenomicRegion> > &boundaries,
    vector<vector<double> > &bins,
    vector< vector<double> > &dead_scales);

void
BinReadsCorrectDeadZones(
    const vector<vector<SimpleGenomicRegion> > &reads, 
    const vector<vector<SimpleGenomicRegion> > &dead_zones,
    const size_t bin_size,
    const double max_dead_proportion,
    const double desert_size,
    vector<SimpleGenomicRegion> &regions,
    vector<vector<SimpleGenomicRegion> > &boundaries,
    vector<vector<double> > &read_counts);

void
BinReadsCorrectDeadZones(
    const vector<SimpleGenomicRegion> &chrom_regions,
    const vector<vector<SimpleGenomicRegion> > &reads, 
    const vector<vector<SimpleGenomicRegion> > &dead_zones,
    const size_t bin_size,
    const double max_dead_proportion,
    const double desert_size,
    vector<SimpleGenomicRegion> &regions,
    vector<vector<SimpleGenomicRegion> > &boundaries,
    vector<vector<double> > &read_counts,
    vector<vector<double> > &dead_scales);


#endif
