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

#ifndef GENOMIC_REGION_TOOL_HPP
#define GENOMIC_REGION_TOOL_HPP

#include "GenomicRegion.hpp"

#include <vector>
#include <string>

using std::vector;
using std::string;

namespace GenomicRegionTool
{
		// input and output operations
		template <class GenomicRegionType> void
		ReadBEDFile(const string filename,
					vector<GenomicRegionType> &genomic_region_set);

		template <class GenomicRegionType> void
		WriteBEDFile(const string filename,
					 vector<GenomicRegionType> &genomic_region_set);

		// functions to get statistical summary of a segmentation scheme
		template <class GenomicRegionType> void
		summary(const vector<GenomicRegionType> &genomic_region_set);
		template <class GenomicRegionType> size_t
		length(const vector<GenomicRegionType> &genomic_region_set);
		   
		// operations on single list of genomics region
		template <class GenomicRegionType> void
		sort(vector<GenomicRegionType> &genomic_region_set);

		void get_contig(const vector<SimpleGenomicRegion> &genomic_region_set,
						vector<SimpleGenomicRegion> &contig_set);
		void get_contig(const vector<GenomicRegion> &genomic_region_set,
						vector<GenomicRegion> &contig_set);

		// seperation and combination operation
		template <class GenomicRegionType> void
		seperate_chromosomes(const vector<GenomicRegionType> &genomic_region_set,
							 vector< vector<GenomicRegionType> > &genomic_region_set_by_chrom);
		template <class GenomicRegionType> void
		seperate(const vector<GenomicRegionType> &genomic_region_set,
				 const vector<GenomicRegionType> &container_region_set,
				 vector< vector<GenomicRegionType> > &result);
		template <class GenomicRegionType> void
		collapse(const vector< vector<GenomicRegionType> > &genomic_region_set_seperate,
				 vector<GenomicRegionType> &genomic_region_set);

		// set operations
		void intersect(const SimpleGenomicRegion &genomic_region_lhs,
					   const SimpleGenomicRegion &genomic_region_rhs,
					   SimpleGenomicRegion &genomic_region_intersect);		
		void intersect(const GenomicRegion &genomic_region_lhs,
					   const GenomicRegion &genomic_region_rhs,
					   GenomicRegion &genomic_region_intersect);		

		void join(const SimpleGenomicRegion &genomic_region_lhs,
				  const SimpleGenomicRegion &genomic_region_rhs,
				  SimpleGenomicRegion &genomic_region_union);
		void join(const GenomicRegion &genomic_region_lhs,
				  const GenomicRegion &genomic_region_rhs,
				  GenomicRegion &genomic_region_union);

		void intersect(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
					   const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
					   vector<SimpleGenomicRegion> &genomic_region_set_intersection);
		void intersect(const vector<GenomicRegion> &genomic_region_set_lhs,
					   const vector<GenomicRegion> &genomic_region_set_rhs,
					   vector<GenomicRegion> &genomic_region_set_intersection);


		template <class GenomicRegionType> void
		join(const vector<GenomicRegionType> &genomic_region_set_lhs,
			 const vector<GenomicRegionType> &genomic_region_set_rhs,
			 vector<GenomicRegionType> &genomic_region_set_union);

		void difference(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
						const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
						vector<SimpleGenomicRegion> &genomic_region_set_difference);

		// comparison of segmentation schemes
		double entropy(const vector<SimpleGenomicRegion> &genomic_region_set,
					   const SimpleGenomicRegion &chrom);
		double mutual_information(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
								  const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
								  const SimpleGenomicRegion &chrom);
		double score_of_similarity(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
								   const vector<SimpleGenomicRegion> &genomic_region_set_rhs);
		double intersection2union_ratio(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
										const vector<SimpleGenomicRegion> &genomic_region_set_rhs);
		double distance(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
						const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
						const SimpleGenomicRegion &chrom);
}

#endif
