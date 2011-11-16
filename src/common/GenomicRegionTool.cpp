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

#include "GenomicRegionTool.hpp"

#include <algorithm>
#include <string>
#include <limits>
#include <iostream>
#include <fstream>
#include <list>

#include <cmath>

using std::vector;
using std::list;


template <class GenomicRegionType> void
GenomicRegionTool::summary(const vector<GenomicRegionType> &genomic_region_set)
{
		const size_t set_size = genomic_region_set.size();
		const size_t total_length = GenomicRegionTool::length(genomic_region_set);
		
		vector<size_t> segment_width_set(set_size);
		for(size_t i = 0; i < set_size; i++)
				segment_width_set[i] = genomic_region_set[i].get_width();
		
		std::sort(segment_width_set.begin(), segment_width_set.end());

		std::cout << std::endl;
		std::cout << "# of segments = "
				  <<  set_size << std::endl;
		std::cout << "total length = "
				  << total_length << std::endl;
		std::cout << "averge segment length = "
				  << (total_length / set_size) << std::endl;
		std::cout << "minimum segment length = "
				  << segment_width_set[0] << std::endl;
		std::cout << "1st quantile = "
				  << segment_width_set[static_cast<size_t>(set_size * 0.25)] << std::endl;
		std::cout << "median segment length = "
				  << segment_width_set[static_cast<size_t>(set_size * 0.50)] << std::endl;
		std::cout << "3rd quantile = "
				  << segment_width_set[static_cast<size_t>(set_size * 0.75)] << std::endl;
		std::cout << "maximum segment length = "
				  << segment_width_set.back() << std::endl;
		std::cout << std::endl;
}

template void
GenomicRegionTool::summary(const vector<GenomicRegion> &genomic_region_set);

template void
GenomicRegionTool::summary(const vector<SimpleGenomicRegion> &genomic_region_set);


/////////////
template <class GenomicRegionType> size_t 
GenomicRegionTool::length(const vector<GenomicRegionType> &genomic_region_set)
{
		size_t s = 0;

		for(size_t i = 0; i < genomic_region_set.size(); ++i)
				s += genomic_region_set[i].get_width();
		return s;
}

template size_t 
GenomicRegionTool::length(const vector<GenomicRegion> &genomic_region_set);

template size_t 
GenomicRegionTool::length(const vector<SimpleGenomicRegion> &genomic_region_set);

/////////////
template <class GenomicRegionType> inline bool 
is_preceed(const GenomicRegionType &lhs, const GenomicRegionType &rhs)
{
		return (lhs.get_chrom() < rhs.get_chrom()) ||
				(lhs.get_chrom() == rhs.get_chrom()
				 && lhs.get_start() < rhs.get_start()) ||
				(lhs.get_chrom() == rhs.get_chrom() 
				 && lhs.get_start() == rhs.get_start()
				 && lhs.get_end() < rhs.get_end());
}

template <class GenomicRegionType> void 
GenomicRegionTool::sort(vector<GenomicRegionType> &genomic_region_set)
{
		std::sort(genomic_region_set.begin(), genomic_region_set.end());
		return;
}
template void 
GenomicRegionTool::sort(vector<GenomicRegion> &genomic_region_set);

template void 
GenomicRegionTool::sort(vector<SimpleGenomicRegion> &genomic_region_set);


/////////////////////

// genomic_region_set should be sorted before passed as arguement
void
GenomicRegionTool::get_contig(const vector<SimpleGenomicRegion> &genomic_region_set,
							  vector<SimpleGenomicRegion> &contig_set)
{
		contig_set.clear();
		
		SimpleGenomicRegion tmp_segment;
		vector<SimpleGenomicRegion>::const_reverse_iterator ri = genomic_region_set.rbegin();
		contig_set.push_back(*ri);
		ri++;
		while(ri != genomic_region_set.rend())
		{
				GenomicRegionTool::join(contig_set.back(),
										*ri,
										tmp_segment);
				if(tmp_segment.get_width() != 0)
						contig_set.back() = tmp_segment;
				else
						contig_set.push_back(*ri);
				ri++;
		}

		std::reverse(contig_set.begin(), contig_set.end());
		return;
}

void
GenomicRegionTool::get_contig(const vector<GenomicRegion> &genomic_region_set,
							  vector<GenomicRegion> &contig_set)
{
/// Here we need to take care of the issue of two strands
		GenomicRegion tmp_segment;
		size_t pos_p_strand = 0;
		size_t pos_n_strand = 0;

		vector<GenomicRegion>::const_reverse_iterator ri = genomic_region_set.rbegin();
		contig_set.push_back(*ri);
		ri++;
		while(ri != genomic_region_set.rend())
		{
				if(ri->get_strand() == '+')
				{
						GenomicRegionTool::join(contig_set[pos_p_strand],
												*ri,
												tmp_segment);
						if(tmp_segment.get_width() != 0)
								contig_set[pos_p_strand] = tmp_segment;
						else
						{
								pos_p_strand = contig_set.size();
								contig_set.push_back(*ri);
						}
				}
				else
				{
						GenomicRegionTool::join(contig_set[pos_n_strand],
												*ri,
												tmp_segment);
						if(tmp_segment.get_width() != 0)
								contig_set[pos_n_strand] = tmp_segment;
						else
						{
								pos_n_strand = contig_set.size();
								contig_set.push_back(*ri);
						}

				}
				ri++;
		}

		std::reverse(contig_set.begin(), contig_set.end());
		return;
}


// input and output functions
// seperated from the original GenomicRegion.hpp
static bool 
is_header_line(const string& line) 
{
		static const char *browser_label = "browser";
		static const size_t browser_label_len = 7;
		for (size_t i = 0; i < browser_label_len; ++i)
				if (line[i] != browser_label[i])
						return false;
		return true;
}


static bool 
is_track_line(const char *line) 
{
		static const char *track_label = "track";
		static const size_t track_label_len = 5;
		for (size_t i = 0; i < track_label_len; ++i)
				if (line[i] != track_label[i])
						return false;
		return true;
}

template <class GenomicRegionType> void
GenomicRegionTool::ReadBEDFile(const string filename,
							   vector<GenomicRegionType> &genomic_region_set)
{
		try
		{
				static const size_t buffer_size = 1000; // Magic
				
				// open and check the file
				std::ifstream in(filename.c_str());
				while (!in.eof()) 
				{
						char buffer[buffer_size];
						in.getline(buffer, buffer_size);
						if (in.gcount() == buffer_size - 1)
								throw;
						if(!is_header_line(buffer) && !is_track_line(buffer))
								genomic_region_set.push_back(GenomicRegionType(buffer));
						in.peek();
				}
				in.close();
		}
		catch(...)
		{
				throw;
		}
		return;

}
template void
GenomicRegionTool::ReadBEDFile<GenomicRegion>(const string filename,
											  vector<GenomicRegion> &genomic_region_set);

template void
GenomicRegionTool::ReadBEDFile<SimpleGenomicRegion>(const string filename,
													vector<SimpleGenomicRegion> &genomic_region_set);

/////////////
template <class GenomicRegionType> void 
GenomicRegionTool::WriteBEDFile(const string filename,
								vector<GenomicRegionType> &genomic_region_set)
{
		try
		{
				std::ofstream out_file(filename.c_str());
				std::copy(genomic_region_set.begin(),
						  genomic_region_set.end(),
						  std::ostream_iterator<GenomicRegionType>(out_file, "\n"));
				out_file.close();
		}
		catch(...)
		{
				throw;
		}
		return;
}

template void 
GenomicRegionTool::WriteBEDFile<SimpleGenomicRegion>(
		const string filename,
		vector<SimpleGenomicRegion> &genomic_region_set);

template void 
GenomicRegionTool::WriteBEDFile<GenomicRegion>(
		const string filename,
		vector<GenomicRegion> &genomic_region_set);

// seperation and combination operation
template <class GenomicRegionType> void 
GenomicRegionTool::seperate_chromosomes(const vector<GenomicRegionType> &genomic_region_set,
										vector< vector<GenomicRegionType> > &genomic_region_set_by_chrom)
{
		genomic_region_set_by_chrom.clear();
		vector<GenomicRegionType> tmp_region_set(genomic_region_set);
		GenomicRegionTool::sort(tmp_region_set);
		
		typename vector<GenomicRegionType>::const_iterator first = tmp_region_set.begin();
		typename vector<GenomicRegionType>::const_iterator last = first;
		while(first != tmp_region_set.end())
		{
				while(last != tmp_region_set.end() &&
					  last->get_chrom() == first->get_chrom())
						last++;
				genomic_region_set_by_chrom.push_back(
						vector<GenomicRegionType>(first, last));

				first = last;
				last = first;
		}
		return;
}
template void 
GenomicRegionTool::seperate_chromosomes<GenomicRegion>(
		const vector<GenomicRegion> &genomic_region_set,
		vector< vector<GenomicRegion> > &genomic_region_set_by_chrom);

template void 
GenomicRegionTool::seperate_chromosomes<SimpleGenomicRegion>(
		const vector<SimpleGenomicRegion> &genomic_region_set,
		vector< vector<SimpleGenomicRegion> > &genomic_region_set_by_chrom);

/////////////////
template <class GenomicRegionType> void
GenomicRegionTool::seperate(const vector<GenomicRegionType> &genomic_region_set,
							const vector<GenomicRegionType> &container_region_set,
							vector< vector<GenomicRegionType> > &result)
{
		result.clear();

		for(size_t i = 0; i < container_region_set.size(); i++)
				result.push_back(vector<GenomicRegionType>());
		
		size_t index = 0;
		size_t container_index = 0;
		while(index < genomic_region_set.size() 
			  && container_index < container_region_set.size())
		{
				

				if(container_region_set[container_index].contains(
						   genomic_region_set[index]))
				{
						result[container_index].push_back(
								genomic_region_set[index]);
				}
				
				if(container_region_set[container_index]
				   < genomic_region_set[index])
						container_index++;
				else
						index++;
		}
		return;

}
template void
GenomicRegionTool::seperate<GenomicRegion>(
		const vector<GenomicRegion> &genomic_region_set,
		const vector<GenomicRegion> &container_region_set,
		vector< vector<GenomicRegion> > &result);

template void
GenomicRegionTool::seperate<SimpleGenomicRegion>(
		const vector<SimpleGenomicRegion> &genomic_region_set,
		const vector<SimpleGenomicRegion> &container_region_set,
		vector< vector<SimpleGenomicRegion> > &result);


/////////////////
template <class GenomicRegionType> void 
GenomicRegionTool::collapse(const vector< vector<GenomicRegionType> > &genomic_region_set_seperate,
							vector<GenomicRegionType> &genomic_region_set)
{
		genomic_region_set.clear();
		for(size_t i = 0; i < genomic_region_set_seperate.size(); i++)
				genomic_region_set.insert(genomic_region_set.end(),
										  genomic_region_set_seperate[i].begin(),
										  genomic_region_set_seperate[i].end());
		GenomicRegionTool::sort(genomic_region_set);
		return;
}
template void 
GenomicRegionTool::collapse<GenomicRegion>(
		const vector< vector<GenomicRegion> > &genomic_region_set_seperate,
		vector<GenomicRegion> &genomic_region_set);

template void 
GenomicRegionTool::collapse<SimpleGenomicRegion>(
		const vector< vector<SimpleGenomicRegion> > &genomic_region_set_seperate,
		vector<SimpleGenomicRegion> &genomic_region_set);


//////////////// set operations
void 
GenomicRegionTool::intersect(const SimpleGenomicRegion &genomic_region_lhs,
							 const SimpleGenomicRegion &genomic_region_rhs,
							 SimpleGenomicRegion &genomic_region_intersect)
{
		if(genomic_region_lhs.get_chrom() == genomic_region_rhs.get_chrom())
		{
				size_t start = std::max(genomic_region_lhs.get_start(),
										genomic_region_rhs.get_start());
				size_t end = std::min(genomic_region_lhs.get_end(),
									  genomic_region_rhs.get_end());
				if(start < end)
						genomic_region_intersect = 
								SimpleGenomicRegion(genomic_region_lhs.get_chrom(),
													start,
													end);
				else
						genomic_region_intersect = 
								SimpleGenomicRegion();
		}
		else
		{
				genomic_region_intersect = 
						SimpleGenomicRegion();
		}
		return;
}

void 
GenomicRegionTool::intersect(const GenomicRegion &genomic_region_lhs,
							 const GenomicRegion &genomic_region_rhs,
							 GenomicRegion &genomic_region_intersect)
{
		const size_t max_name_len = 128;
		
		if(genomic_region_lhs.get_chrom() == genomic_region_rhs.get_chrom()
		   && genomic_region_rhs.get_strand() == genomic_region_rhs.get_strand())
		{
				size_t start = std::max(genomic_region_lhs.get_start(),
										genomic_region_rhs.get_start());
				size_t end = std::min(genomic_region_lhs.get_end(),
									  genomic_region_rhs.get_end());
				if(start < end)
				{
						float score = (genomic_region_lhs.get_score()
									   + genomic_region_rhs.get_score()) / 2;
						string name("");
						name += genomic_region_lhs.get_name().substr(0, max_name_len / 2);
						name += genomic_region_rhs.get_name().substr(0, max_name_len / 2);
						
						genomic_region_intersect = 
								GenomicRegion(genomic_region_lhs.get_chrom(),
											  start,
											  end,
											  name,
											  score,
											  genomic_region_lhs.get_strand());
				}
				else
						genomic_region_intersect = 
								GenomicRegion();
		}
		else
		{
				genomic_region_intersect = 
						GenomicRegion();
		}
		return;
}

///////////////////

void 
GenomicRegionTool::join(const SimpleGenomicRegion &genomic_region_lhs,
						const SimpleGenomicRegion &genomic_region_rhs,
						SimpleGenomicRegion &genomic_region_union)
{
		if(genomic_region_lhs.get_chrom() == genomic_region_rhs.get_chrom())
		{
				size_t start = std::min(genomic_region_lhs.get_start(),
										genomic_region_rhs.get_start());
				size_t end = std::max(genomic_region_lhs.get_end(),
									  genomic_region_rhs.get_end());
				if(end  - start <=
				   genomic_region_lhs.get_width() + genomic_region_rhs.get_width())
						genomic_region_union = 
								SimpleGenomicRegion(genomic_region_lhs.get_chrom(),
													start,
													end);
				else
						genomic_region_union = 
								SimpleGenomicRegion();
		}
		else
		{
				genomic_region_union = 
						SimpleGenomicRegion();
		}
		return;
}

void 
GenomicRegionTool::join(const GenomicRegion &genomic_region_lhs,
						const GenomicRegion &genomic_region_rhs,
						GenomicRegion &genomic_region_union)
{
		const size_t max_name_len = 128;

		if(genomic_region_lhs.get_chrom() == genomic_region_rhs.get_chrom()
		   && genomic_region_rhs.get_strand() == genomic_region_rhs.get_strand())
		{
				size_t start = std::min(genomic_region_lhs.get_start(),
										genomic_region_rhs.get_start());
				size_t end = std::max(genomic_region_lhs.get_end(),
									  genomic_region_rhs.get_end());
				if(end  - start <=
				   genomic_region_lhs.get_width() + genomic_region_rhs.get_width())

				{
						float score = (genomic_region_lhs.get_score()
									   + genomic_region_rhs.get_score()) / 2;
						string name("");
						name += genomic_region_lhs.get_name().substr(0, max_name_len / 2);
						name += genomic_region_rhs.get_name().substr(0, max_name_len / 2);
						
						genomic_region_union = 
								GenomicRegion(genomic_region_lhs.get_chrom(),
											  start,
											  end,
											  name,
											  score,
											  genomic_region_lhs.get_strand());
				}
				else
						genomic_region_union = 
								GenomicRegion();
		}
		else
		{
				genomic_region_union = 
						GenomicRegion();
		}
		return;
}



// CAUTION:
// THIS VERSION ONLY WORKS FOR UNORVERLAPPING SEGMENTS SERIES
// OR ALL SEGMENTS HAVE EQUAL LENGTH
void 
GenomicRegionTool::intersect(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
							 const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
							 vector<SimpleGenomicRegion> &genomic_region_set_intersection)
{
		genomic_region_set_intersection.clear();

		const size_t lhs_size = genomic_region_set_lhs.size();
		const size_t rhs_size = genomic_region_set_rhs.size();

		size_t lhs_counter = 0, rhs_counter = 0;
		SimpleGenomicRegion tmp_segment;

		while(lhs_counter < lhs_size && rhs_counter < rhs_size)
		{
				GenomicRegionTool::intersect(genomic_region_set_lhs[lhs_counter],
											 genomic_region_set_rhs[rhs_counter],
											 tmp_segment);
				if(tmp_segment.get_width() > 0)
						genomic_region_set_intersection.push_back(tmp_segment);
				if(genomic_region_set_lhs[lhs_counter]
				   <  genomic_region_set_rhs[rhs_counter])
						lhs_counter++;
				else
						rhs_counter++;
		}
		return;
}

void 
GenomicRegionTool::intersect(const vector<GenomicRegion> &genomic_region_set_lhs,
							 const vector<GenomicRegion> &genomic_region_set_rhs,
							 vector<GenomicRegion> &genomic_region_set_intersection)
{
		///// not well defined
		return;
}

//////////////////
template <class GenomicRegionType> void
GenomicRegionTool::join(const vector<GenomicRegionType> &genomic_region_set_lhs,
						const vector<GenomicRegionType> &genomic_region_set_rhs,
						vector<GenomicRegionType> &genomic_region_set_union)
{
		vector<GenomicRegionType> tmp_region_set;
		typename vector<GenomicRegionType>::const_iterator lhs_iterator =
				genomic_region_set_lhs.begin();
		typename vector<GenomicRegionType>::const_iterator rhs_iterator =
				genomic_region_set_rhs.begin();

		while(lhs_iterator != genomic_region_set_lhs.end()
			  && rhs_iterator != genomic_region_set_rhs.end())
		{
				if(*lhs_iterator < *rhs_iterator)
				{
						tmp_region_set.push_back(*lhs_iterator);
						lhs_iterator++;
				}
				else
				{
						tmp_region_set.push_back(*rhs_iterator);
						rhs_iterator++;
				}
		}

		while(lhs_iterator != genomic_region_set_lhs.end())
		{
				tmp_region_set.push_back(*lhs_iterator);
				lhs_iterator++;
		}
		while(rhs_iterator != genomic_region_set_rhs.end())
		{
						tmp_region_set.push_back(*rhs_iterator);
						rhs_iterator++;
		}
		
		GenomicRegionTool::get_contig(tmp_region_set,
									  genomic_region_set_union);
		return;
}

template void
GenomicRegionTool::join<SimpleGenomicRegion>(
		const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
		const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
		vector<SimpleGenomicRegion> &genomic_region_set_union);

template void
GenomicRegionTool::join<GenomicRegion>(
		const vector<GenomicRegion> &genomic_region_set_lhs,
		const vector<GenomicRegion> &genomic_region_set_rhs,
		vector<GenomicRegion> &genomic_region_set_union);


///////////////////

void
GenomicRegionTool::difference(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
		   const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
		   vector<SimpleGenomicRegion> &genomic_region_set_difference)
{
		genomic_region_set_difference.clear();
		
		const size_t max_chrom_length = std::numeric_limits<size_t>::max();

		// get the complement set of right-hand-side genomic_region_set
		vector<SimpleGenomicRegion> complement_set;
		vector<SimpleGenomicRegion>::const_iterator
				lagging_iterator = genomic_region_set_rhs.begin();
		vector<SimpleGenomicRegion>::const_iterator
				leading_iterator = 	genomic_region_set_rhs.begin();

		SimpleGenomicRegion tmp_segment;
		
		tmp_segment.set_chrom(leading_iterator->get_chrom());
		tmp_segment.set_start(0);
		tmp_segment.set_end(leading_iterator->get_start());
		complement_set.push_back(tmp_segment);
		leading_iterator++;
		while(leading_iterator != genomic_region_set_rhs.end())
		{
				if(leading_iterator->get_chrom() == lagging_iterator->get_chrom())
				{
						SimpleGenomicRegion tmp_segment;
						tmp_segment.set_chrom(leading_iterator->get_chrom());
						tmp_segment.set_start(lagging_iterator->get_end());
						tmp_segment.set_end(leading_iterator->get_start());
						complement_set.push_back(tmp_segment);
				}
				else
				{
						SimpleGenomicRegion tmp_segment;
						tmp_segment.set_chrom(lagging_iterator->get_chrom());
						tmp_segment.set_start(lagging_iterator->get_end());
						tmp_segment.set_end(max_chrom_length);
						complement_set.push_back(tmp_segment);

						tmp_segment.set_chrom(leading_iterator->get_chrom());
						tmp_segment.set_start(0);
						tmp_segment.set_end(leading_iterator->get_start());
						complement_set.push_back(tmp_segment);
				}
				lagging_iterator = leading_iterator;
				leading_iterator++;
		}
		complement_set.push_back(SimpleGenomicRegion(lagging_iterator->get_chrom(),
													 lagging_iterator->get_end(),
													 max_chrom_length));
				
		GenomicRegionTool::intersect(genomic_region_set_lhs,
									 complement_set,
									 genomic_region_set_difference);
		return;
}

// add an implementation of score_of_similarity
// but it seems not reasonable in some cases: 
// for example, if two identical set of segmentation schemes,
// it gives much more score to those who identify more segments as foreground
double 
GenomicRegionTool::score_of_similarity(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
									   const vector<SimpleGenomicRegion> &genomic_region_set_rhs)
{
		double episilon = 1e-6;
		
		std::vector<SimpleGenomicRegion> tmp_region_set;
		GenomicRegionTool::intersect(genomic_region_set_lhs, genomic_region_set_rhs, tmp_region_set);
		size_t inter_length = GenomicRegionTool::length(tmp_region_set);
		
		GenomicRegionTool::difference(genomic_region_set_lhs, genomic_region_set_rhs, tmp_region_set);
		size_t diff_len = GenomicRegionTool::length(tmp_region_set);

		GenomicRegionTool::difference(genomic_region_set_rhs, genomic_region_set_lhs, tmp_region_set);
		diff_len += GenomicRegionTool::length(tmp_region_set);
		
		return inter_length / (diff_len + episilon);
}


// compute the entropy of an segmentation
// assume single chromosome
inline double
log_base(const double v, const double base = 2.0)
{
		return std::log(v) / std::log(base);
}

double 
GenomicRegionTool::entropy(const vector<SimpleGenomicRegion> &genomic_region_set,
			   const SimpleGenomicRegion &chrom)
{
		const double chrom_len = static_cast<double>(chrom.get_width());
		double h = 0; 			// entropy
		double p = 0;			// probability
		
		size_t last_end = 0;
		for(size_t i = 0; i < genomic_region_set.size(); i++)
		{
				p = (genomic_region_set[i].get_start() - last_end) / chrom_len;
				h += - p * log_base(p);
				
				p = genomic_region_set[i].get_width() / chrom_len;
				h += - p * log_base(p);
				
				last_end = genomic_region_set[i].get_end();
		}
 		p = (chrom.get_end() - last_end) / chrom_len;
 		h += (chrom.get_end() > last_end) ? (- p * log_base(p)) : 0;

		return h;
}

// compute the mutual information between segmentation schemes
// we assume each segmentation only deals with one chromorome
double 
GenomicRegionTool::mutual_information(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
									  const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
									  const SimpleGenomicRegion &chrom)
{
		vector<SimpleGenomicRegion> full_set_lhs, full_set_rhs;
		vector<SimpleGenomicRegion> chrom_set;
		chrom_set.push_back(chrom);
		GenomicRegionTool::difference(chrom_set,
									  genomic_region_set_lhs,
									  full_set_lhs);
		full_set_lhs.insert(full_set_lhs.end(),
							genomic_region_set_lhs.begin(),
							genomic_region_set_lhs.end());
		GenomicRegionTool::sort(full_set_lhs);
		
		GenomicRegionTool::difference(chrom_set,
									  genomic_region_set_rhs,
									  full_set_rhs);
		full_set_rhs.insert(full_set_rhs.end(),
							genomic_region_set_rhs.begin(),
							genomic_region_set_rhs.end());
		GenomicRegionTool::sort(full_set_rhs);

		const double chrom_len = static_cast<double>(chrom.get_width());
		double mutual_info = 0; 			// mutual information
		SimpleGenomicRegion tmp_intersection;
		size_t tmp_intersection_len;

		const size_t lhs_size = full_set_lhs.size();
		const size_t rhs_size = full_set_rhs.size();
		
		size_t lhs_counter = 0;
		size_t rhs_counter = 0;
		while(lhs_counter < lhs_size && rhs_counter < rhs_size)
		{
				GenomicRegionTool::intersect(full_set_lhs[lhs_counter],
											 full_set_rhs[rhs_counter],
											 tmp_intersection);

				if((tmp_intersection_len = tmp_intersection.get_width()) > 0)
						mutual_info += (tmp_intersection_len / chrom_len)
								* log_base(tmp_intersection_len 
										   * chrom_len 
										   / full_set_lhs[lhs_counter].get_width()
										   / full_set_rhs[rhs_counter].get_width());
						
				if(full_set_lhs[lhs_counter].get_end()
				   < full_set_rhs[rhs_counter].get_end())
						lhs_counter++;
				else
						rhs_counter++;
		}
		
		return mutual_info;
}

double 
GenomicRegionTool::distance(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
							const vector<SimpleGenomicRegion> &genomic_region_set_rhs,
							const SimpleGenomicRegion &chrom)
{
// 		vector<SimpleGenomicRegion> tmp_region_set;
// 		GenomicRegionTool::intersect(genomic_region_set_lhs,
// 									 genomic_region_set_rhs,
// 									 tmp_region_set);
// 		size_t inter_length = GenomicRegionTool::length(tmp_region_set);
// 		GenomicRegionTool::difference(genomic_region_set_lhs,
// 									  genomic_region_set_rhs,
// 									  tmp_region_set);
// 		size_t diff_length = GenomicRegionTool::length(tmp_region_set);
		
		double mutual_info =
				GenomicRegionTool::mutual_information(genomic_region_set_lhs,
													  genomic_region_set_rhs,
													  chrom);
		double entro =  std::max(GenomicRegionTool::entropy(genomic_region_set_lhs,
															chrom),
								 GenomicRegionTool::entropy(genomic_region_set_rhs,
															chrom));
		return 1 - mutual_info / entro;
}



double
GenomicRegionTool:: intersection2union_ratio(const vector<SimpleGenomicRegion> &genomic_region_set_lhs,
											 const vector<SimpleGenomicRegion> &genomic_region_set_rhs)
{
		std::vector<SimpleGenomicRegion> tmp_region_set;
		GenomicRegionTool::intersect(genomic_region_set_lhs, 
									 genomic_region_set_rhs, 
									 tmp_region_set);
		size_t intersect_set_len = GenomicRegionTool::length(tmp_region_set);

		GenomicRegionTool::join(genomic_region_set_lhs,
								genomic_region_set_rhs,
								tmp_region_set);
		size_t union_set_len = GenomicRegionTool::length(tmp_region_set);

		return static_cast<double>(intersect_set_len) / union_set_len;
}



