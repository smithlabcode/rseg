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

#ifndef SELECT_BIN_SIZE
#define SELECT_BIN_SIZE

#include <vector>

#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

size_t
select_bin_size(const size_t n_reads, const size_t genome_size,
		const double alpha, const double theta);

size_t
select_bin_size_naive(const std::vector<SimpleGenomicRegion> &regions,
                      const std::vector<std::vector<SimpleGenomicRegion> > &reads,
                      const std::vector<std::vector<SimpleGenomicRegion> > &deads);

size_t
select_bin_size_waterman(const std::vector<SimpleGenomicRegion> &regions,
                         const std::vector<std::vector<SimpleGenomicRegion> > &reads,
                         const std::vector<std::vector<SimpleGenomicRegion> > &deads,
                         const bool smooth = true);

size_t
select_bin_size_hideaki(const std::vector<SimpleGenomicRegion> &regions,
                        const std::vector<std::vector<SimpleGenomicRegion> > &reads,
                        const std::vector<std::vector<SimpleGenomicRegion> > &deads,
                        const bool smooth = true);
size_t 
select_bin_size_hideaki_emp(const std::vector<SimpleGenomicRegion> &regions,
                            const std::vector<std::vector<SimpleGenomicRegion> > &reads,
                            const std::vector<std::vector<SimpleGenomicRegion> > &deads,
                            const double max_dead_proportion);
#endif
