/* 
 * Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
 *
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

#include "ReadCounts.hpp"

#include "RNG.hpp"
#include <cmath>

using std::vector;
using std::string;
using std::min;
using std::max;

using std::cerr;
using std::endl;

// static void
// BinReads(const vector<SimpleGenomicRegion> &reads,
// 	 const size_t region_start, 
// 	 const size_t region_end,
// 	 const size_t bin_size, 
// 	 vector<double> &bins) 
// {
//   bins.clear();
//   size_t read_idx = 0;
//   for (size_t i = region_start; i <= region_end - bin_size; i += bin_size) {
//     size_t counts = 0;
//     while (read_idx < reads.size() && reads[read_idx].get_start() < i + bin_size) {
//       if (reads[read_idx].get_start() >= i)
// 	++counts;
//       ++read_idx;
//     }
//     bins.push_back(counts);
//   }
// }


// void
// BinReads(const vector<vector<SimpleGenomicRegion> > &reads,
// 	 const vector<SimpleGenomicRegion> &regions,
// 	 const size_t bin_size,
// 	 vector<vector<double> > &bins) {
//   assert(reads.size() == regions.size());
//   bins.resize(reads.size());
//   for (size_t i = 0; i < regions.size(); ++i) {
//     BinReads(reads[i], regions[i].get_start(), regions[i].get_end(),
// 	     bin_size, bins[i]);
//   }
// }



// static void
// BinReadsExtendOverDeadZones(const vector<SimpleGenomicRegion> &reads, 
// 			    const vector<SimpleGenomicRegion> &dead_zones,
// 			    const string chrom_name,
// 			    const size_t region_start,
// 			    const size_t region_end,
// 			    const size_t bin_size,
// 			    vector<SimpleGenomicRegion> &boundaries,
// 			    vector<double> &bins) {
  
//   bins.clear();
//   boundaries.clear();
  
//   vector<SimpleGenomicRegion>::const_iterator ditr(dead_zones.begin());
//   const vector<SimpleGenomicRegion>::const_iterator dlim(dead_zones.end());
  
//   size_t read_idx = 0;
//   size_t i = region_start;
//   while (i < region_end) {
//     size_t end = i;
//     size_t goods = 0;
    
//     if (ditr != dlim && end + bin_size > ditr->get_start()) {
//       while (goods < bin_size) {
// 	if (ditr < dlim) {
	  
// 	  const size_t curr_dead_end = ditr->get_end();
// 	  const size_t curr_dead_start = ditr->get_start();
	  
// 	  if (end >= curr_dead_end)
// 	    ++ditr;
// 	  else {
// 	    if (end < curr_dead_start) {
// 	      const size_t new_goods = min(bin_size - goods, curr_dead_start - end);
// 	      goods += new_goods;
// 	      end += new_goods;
// 	    }
// 	    else end = curr_dead_end;
// 	  }
// 	}
// 	else {
// 	  end += bin_size - goods;
// 	  goods = bin_size;
// 	}
//       }
//     }
//     else end += bin_size;
    
//     size_t counts = 0;
//     while (read_idx < reads.size() && reads[read_idx].get_start() < end) {
//       if (reads[read_idx].get_start() >= i)
// 	counts++;
//       read_idx++;
//     }
//     bins.push_back(counts);
//     boundaries.push_back(SimpleGenomicRegion(chrom_name, i, end));
//     i = end;
//   }
// }

// void
// BinReadsExtendOverDeadZones(const vector<vector<SimpleGenomicRegion> > &reads, 
// 			    const vector<vector<SimpleGenomicRegion> > &dead_zones,
// 			    const vector<SimpleGenomicRegion> &regions,
// 			    const size_t bin_size,
// 			    vector<vector<SimpleGenomicRegion> > &boundaries,
// 			    vector<vector<double> > &bins) {
//   assert(reads.size() == regions.size());
//   assert(dead_zones.size() == regions.size());
//   bins.resize(reads.size());
//   boundaries.resize(reads.size());
//   for (size_t i = 0; i < regions.size(); ++i) {
//     BinReadsExtendOverDeadZones(reads[i], 
// 				dead_zones[i], regions[i].get_chrom(), regions[i].get_start(), 
// 				regions[i].get_end(),bin_size, boundaries[i], bins[i]);
//   }
// }

void
BinDeads(const vector<SimpleGenomicRegion> &deads,
	 const size_t region_start, const size_t region_end,
	 const size_t bin_size, vector<double> &bins) {
  bins.clear();
  size_t dead_idx = 0;
  for (size_t i = region_start; i < region_end; i += bin_size) {
    size_t counts = 0;
    while (dead_idx < deads.size() && deads[dead_idx].get_end() < i + bin_size) {
      /* Don't need the condition below because the start is always
	 within the first bin, since that's how the regions were
	 determined, and the bins start at the start of the first
	 region. */
      // if (deads[dead_idx].get_start() >= i)   ?????
      counts += deads[dead_idx].get_end() - max(deads[dead_idx].get_start(), i);
      ++dead_idx;
    }
    if (dead_idx < deads.size() && deads[dead_idx].get_start() < i + bin_size)
      counts += ((i + bin_size) - max(deads[dead_idx].get_start(), i));
    assert(counts <= bin_size);
    bins.push_back(counts);
  }
}

// static void
// BinReadsCorrectDeadZones(const vector<SimpleGenomicRegion> &reads, 
// 			 const vector<SimpleGenomicRegion> &dead_zones,
// 			 const string chrom_name,
// 			 const size_t region_start,
// 			 const size_t region_end,
// 			 const size_t bin_size,
// 			 const double max_dead_proportion,
// 			 vector<SimpleGenomicRegion> &boundaries,
// 			 vector<double> &bins) 
// {
//   vector<double> read_bins;
//   BinReads(reads, region_start, region_end, bin_size, read_bins);
  
//   vector<double> dead_bins;
//   BinDeads(dead_zones, region_start, region_end, bin_size, dead_bins);
  


//   //     vector<double> corr;
//   //     for (size_t i = 0; i < dead_bins.size(); ++i)
//   //         corr.push_back(bin_size/(bin_size - dead_bins[i]));

//   // (SQ) !!! it is possible that  bin_size == dead_bins 
//   vector<double> corr;
//   for (size_t i = 0; i < dead_bins.size(); ++i)
//     if (bin_size != dead_bins[i])
//       corr.push_back(bin_size/(bin_size - dead_bins[i]));
//     else
//       corr.push_back(1.0); 
    

//   Runif rng(time(0) + getpid());
//   boundaries.clear();
//   bins.clear();
//   for (size_t i = 0; i < corr.size(); ++i) {
//     if (dead_bins[i] < max_dead_proportion*bin_size)
//       {
// 	boundaries.push_back(SimpleGenomicRegion(chrom_name,
// 						 region_start + i*bin_size,
// 						 region_start + (i + 1)*bin_size));
// 	const double corrected = read_bins[i]*corr[i];
// 	const double floor_count = std::floor(corrected);
// 	const double frac_part = corrected - floor_count;
// 	bins.push_back((rng.runif(0.0, 1.1) < frac_part) ? 
// 		       std::ceil(corrected) : floor_count);
//       }
//   }
// }

// void
// BinReadsCorrectDeadZones(const vector<vector<SimpleGenomicRegion> > &reads, 
// 			 const vector<vector<SimpleGenomicRegion> > &dead_zones,
// 			 const vector<SimpleGenomicRegion> &regions,
// 			 const size_t bin_size,
// 			 const double max_dead_proportion,
// 			 vector<vector<SimpleGenomicRegion> > &boundaries,
// 			 vector<vector<double> > &bins) 
// {
//   assert(reads.size() == regions.size() && dead_zones.size() == reads.size());
//   bins.resize(reads.size());
//   boundaries.resize(reads.size());
//   for (size_t i = 0; i < regions.size(); ++i)
//     BinReadsCorrectDeadZones(reads[i], dead_zones[i], regions[i].get_chrom(), 
// 			     regions[i].get_start(), regions[i].get_end(),
// 			     bin_size, max_dead_proportion, boundaries[i], bins[i]);
// }

static void
BinReadsCorrectDeadZones(const vector<SimpleGenomicRegion> &dead_zones,
			 const size_t start, const size_t end,
			 const size_t bin_size, vector<double> &dead_scales) {
  vector<double> dead_bins;
  BinDeads(dead_zones, start, end, bin_size, dead_bins);
  dead_scales.clear();
  for (size_t i = 0; i < dead_bins.size(); ++i)
    dead_scales.push_back(1 - dead_bins[i]/bin_size);
}

void
BinReadsCorrectDeadZones(const vector<vector<SimpleGenomicRegion> > &dead_zones,
			 const vector<SimpleGenomicRegion> &regions,
			 const size_t bin_size,
			 vector< vector<double> > &dead_scales) {
  dead_scales.resize(regions.size());
  for (size_t i = 0; i < regions.size(); ++i)
    BinReadsCorrectDeadZones(dead_zones[i], regions[i].get_start(), regions[i].get_end(),
			     bin_size, dead_scales[i]);
}

// ///// 
// static void
// BinReadsCorrectDeadZones(const vector<SimpleGenomicRegion> &reads, 
// 			 const vector<SimpleGenomicRegion> &dead_zones,
// 			 const string chrom_name,
// 			 const size_t region_start,
// 			 const size_t region_end,
// 			 const size_t bin_size,
// 			 const double max_dead_proportion,
//                          const double desert_size,
// 			 vector<vector<SimpleGenomicRegion> > &boundaries,
// 			 vector<vector<double> > &bins) { 
//   vector<double> read_bins;
//   BinReads(reads, region_start, region_end, bin_size, read_bins);
    
//   vector<double> dead_bins;
//   BinDeads(dead_zones, region_start, region_end, bin_size, dead_bins);
  
//   Runif rng(time(0) + getpid());

//   boundaries.resize(1);
//   bins.resize(1);
//   size_t gap_bin_num = 0;
    
//   for (size_t i = 0; i < dead_bins.size(); ++i)
//     if (dead_bins[i] < max_dead_proportion*bin_size) {
//       if (gap_bin_num * bin_size > desert_size) {
// 	boundaries.push_back(vector<SimpleGenomicRegion>());
// 	bins.push_back(vector<double>());
//       }
      
//       gap_bin_num = 0;
      
//       const double corr = bin_size / (bin_size - dead_bins[i]);
//       const double corrected = read_bins[i] * corr;
//       const double floor_count = std::floor(corrected);
//       const double frac_part = corrected - floor_count;
//       bins.back().push_back((rng.runif(0.0, 1.1) < frac_part) ? 
// 			    std::ceil(corrected) : floor_count);
      
//       boundaries.back().push_back(
// 				  SimpleGenomicRegion(chrom_name,
// 						      region_start + i*bin_size,
// 						      region_start + (i + 1)*bin_size));
//     }
//     else
//       ++gap_bin_num;
// }

// void
// BinReadsCorrectDeadZones(const vector<vector<SimpleGenomicRegion> > &reads, 
// 			 const vector<vector<SimpleGenomicRegion> > &dead_zones,
// 			 const size_t bin_size,
// 			 const double max_dead_proportion,
// 			 const double desert_size,
// 			 vector<SimpleGenomicRegion> &regions,
// 			 vector<vector<SimpleGenomicRegion> > &boundaries,
// 			 vector<vector<double> > &read_counts) {
//   assert(reads.size() == dead_zones.size());

//   regions.clear();
//   boundaries.clear();
//   read_counts.clear();
    
//   for (size_t i = 0; i < reads.size(); ++i) {
//     vector<vector<SimpleGenomicRegion> >  tmp_bins;
//     vector<vector<double> > tmp_read_counts;
//     vector<vector<double> > tmp_dead_scales;
        
//     const string chrom = reads[i].front().get_chrom();
//     const size_t start = std::min(reads[i].front().get_start(),
// 				  dead_zones[i].front().get_start());
//     const size_t end = std::max(reads[i].back().get_end(),
// 				dead_zones[i].back().get_end());
        
//     BinReadsCorrectDeadZones(reads[i], dead_zones[i], chrom, start, end,
// 			     bin_size, max_dead_proportion, desert_size,
// 			     tmp_bins, tmp_read_counts);
    
//     for (size_t j = 0; j < tmp_bins.size(); ++j)
//       if (tmp_bins[j].size()) {
// 	regions.push_back(SimpleGenomicRegion(chrom, tmp_bins[j][0].get_start(),
// 					      tmp_bins[j][0].get_end()));
// 	boundaries.push_back(tmp_bins[j]);
// 	read_counts.push_back(tmp_read_counts[j]);
//       }
//   }
// }


// static void
// BinReadsCorrectDeadZones(const vector<SimpleGenomicRegion> &reads, 
// 			 const vector<SimpleGenomicRegion> &dead_zones,
// 			 const string chrom_name,
// 			 const size_t region_start,
// 			 const size_t region_end,
// 			 const size_t bin_size,
// 			 const double max_dead_proportion,
//                          const double desert_size,
// 			 vector<vector<SimpleGenomicRegion> > &boundaries,
// 			 vector<vector<double> > &bins,
//                          vector<vector<double> > &dead_scales) { 
  
//   vector<double> read_bins;
//   BinReads(reads, region_start, region_end, bin_size, read_bins);
  
//   vector<double> dead_bins;
//   BinDeads(dead_zones, region_start, region_end, bin_size, dead_bins);
  
//   boundaries.resize(1);
//   bins.resize(1);
//   dead_scales.resize(1);
//   size_t gap_bin_num = 0;
  
//   for (size_t i = 0; i < dead_bins.size(); ++i) {
//     if (dead_bins[i] < max_dead_proportion*bin_size) {
//       if (gap_bin_num * bin_size > desert_size) {
// 	boundaries.push_back(vector<SimpleGenomicRegion>());
// 	bins.push_back(vector<double>());
// 	dead_scales.push_back(vector<double>());
//       }
      
//       gap_bin_num = 0;
//       boundaries.back().push_back(SimpleGenomicRegion(chrom_name,
// 						      region_start + i*bin_size,
// 						      region_start + (i + 1)*bin_size));
//       bins.back().push_back(read_bins[i]);
//       dead_scales.back().push_back(1 - dead_bins[i]/bin_size);
//     }
//     else ++gap_bin_num;
//   }
// }

// void
// BinReadsCorrectDeadZones(const vector<SimpleGenomicRegion> &chrom_regions,
// 			 const vector<vector<SimpleGenomicRegion> > &reads, 
// 			 const vector<vector<SimpleGenomicRegion> > &dead_zones,
// 			 const size_t bin_size,
// 			 const double max_dead_proportion,
// 			 const double desert_size,
// 			 vector<SimpleGenomicRegion> &regions,
// 			 vector<vector<SimpleGenomicRegion> > &boundaries,
// 			 vector<vector<double> > &read_counts,
// 			 vector<vector<double> > &dead_scales) {
  
//   assert(reads.size() == dead_zones.size());
  
//   regions.clear();
//   boundaries.clear();
//   read_counts.clear();
//   dead_scales.clear();
  
//   for (size_t i = 0; i < chrom_regions.size(); ++i) {
//     vector<vector<SimpleGenomicRegion> >  tmp_bins;
//     vector<vector<double> > tmp_read_counts;
//     vector<vector<double> > tmp_dead_scales;
    
//     const string chrom(chrom_regions[i].get_chrom());
//     const size_t start = chrom_regions[i].get_start();
//     const size_t end = chrom_regions[i].get_end();
        
//     BinReadsCorrectDeadZones(reads[i], dead_zones[i], chrom, start, end,
// 			     bin_size, max_dead_proportion, desert_size,
// 			     tmp_bins, tmp_read_counts, tmp_dead_scales);
    
//     for (size_t j = 0; j < tmp_bins.size(); ++j)
//       if (tmp_bins[j].size()) {
// 	regions.push_back(SimpleGenomicRegion(chrom,
// 					      tmp_bins[j].front().get_start(),
// 					      tmp_bins[j].back().get_end()));
// 	boundaries.push_back(tmp_bins[j]);
// 	read_counts.push_back(tmp_read_counts[j]);
// 	dead_scales.push_back(tmp_dead_scales[j]);
//       }
//   }
// }
