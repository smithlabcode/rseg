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

#include "LoadReadsByRegion.hpp"
#include "SortGenomicRegion.hpp"

#include "smithlab_utils.hpp"

#include <iomanip>

using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::setw;
using std::min;
using std::max;



static void
get_deserts(const vector<SimpleGenomicRegion> &dead_zones, 
			const size_t desert_size,
			vector<SimpleGenomicRegion> &deserts) 
{
    deserts.clear();
    for (vector<SimpleGenomicRegion>::const_iterator i = 
             dead_zones.begin(); i != dead_zones.end(); ++i)
        if (i->get_width() > desert_size)
            deserts.push_back(*i);
}



static void
separate_regions(const vector<SimpleGenomicRegion> &big_regions,
				 const vector<SimpleGenomicRegion> &regions, 
				 vector<vector<SimpleGenomicRegion> > &sep_regions) 
{
    size_t rr_id = 0;
    const size_t n_regions = regions.size();
    const size_t n_big_regions = big_regions.size();
    sep_regions.resize(n_big_regions);
    for (size_t i = 0; i < n_big_regions; ++i) 
    {
        const string current_chrom(big_regions[i].get_chrom());
        const size_t current_start = big_regions[i].get_start();
        const size_t current_end = big_regions[i].get_end();
        while (rr_id < n_regions &&
               (regions[rr_id].get_chrom() < current_chrom ||
                (regions[rr_id].get_chrom() == current_chrom &&
                 regions[rr_id].get_start() < current_start)))
            ++rr_id;
        while (rr_id < n_regions &&
               (regions[rr_id].get_chrom() == current_chrom &&
                regions[rr_id].get_start() < current_end)) 
        {
            sep_regions[i].push_back(regions[rr_id]);
            ++rr_id;
        }
    }
}


// static void
// separate_regions(const vector<SimpleGenomicRegion> &big_regions,
// 				 const vector<SimpleGenomicRegion> &regions, 
// 				 vector<vector<SimpleGenomicRegion> > &sep_regions) 
// {
//     size_t rr_id = 0;
//     const size_t n_regions = regions.size();
//     const size_t n_big_regions = big_regions.size();
//     sep_regions.resize(n_big_regions);
//     for (size_t i = 0; i < n_big_regions; ++i) 
//     {
//         while (rr_id < n_regions && (! big_regions[i].contains(regions[rr_id])))
//             ++rr_id;

//         while (rr_id < n_regions && big_regions[i].contains(regions[rr_id]))
//         {
//             sep_regions[i].push_back(regions[rr_id]);
//             ++rr_id;
//         }
//     }
// }



/* This function inverts a set of desert regions with respect to a set
   of regions.
*/
static void
invert_deserts(const vector<SimpleGenomicRegion> &regions,
			   const vector<SimpleGenomicRegion> &deserts,
			   vector<SimpleGenomicRegion> &inverted) 
{

    // separate the deserts by region, truncating them at the ends
    size_t n_regions = regions.size();
    size_t n_deserts = deserts.size();
    size_t desert_id = 0;
    vector<vector<SimpleGenomicRegion> > sep_deserts(n_regions);
    for (size_t i = 0; i < n_regions; ++i) 
    {
        const string current_chrom(regions[i].get_chrom());
        const size_t current_start = regions[i].get_start();
        const size_t current_end = regions[i].get_end();
        while (desert_id < n_deserts &&
               (deserts[desert_id].get_chrom() < current_chrom ||
                (deserts[desert_id].get_chrom() == current_chrom &&
                 deserts[desert_id].get_end() <= current_start)))
            ++desert_id;
        while (desert_id < n_deserts &&
               (deserts[desert_id].get_chrom() == current_chrom &&
                deserts[desert_id].get_start() < current_end)) 
        {
            SimpleGenomicRegion d(deserts[desert_id]);
            d.set_start(max(d.get_start(), regions[i].get_start()));
            d.set_end(min(d.get_end(), regions[i].get_end()));
            sep_deserts[i].push_back(d);
            ++desert_id;
        }
    }
  
    for (size_t i = 0; i < sep_deserts.size(); ++i) 
    {
        if (!sep_deserts[i].empty()) 
        {
            const string chrom(sep_deserts[i].front().get_chrom());
            if (sep_deserts[i].front().get_start() > regions[i].get_start())
                inverted.push_back(SimpleGenomicRegion(chrom, regions[i].get_start(),
                                                       sep_deserts[i].front().get_start()));
            if (sep_deserts[i].size() > 1)
                for (size_t j = 0; j < sep_deserts[i].size() - 1; ++j)
                    inverted.push_back(SimpleGenomicRegion(chrom, sep_deserts[i][j].get_end(),
                                                           sep_deserts[i][j + 1].get_start()));
            if (sep_deserts[i].back().get_end() < regions[i].get_end())
                inverted.push_back(SimpleGenomicRegion(chrom, sep_deserts[i].back().get_end(),
                                                       regions[i].get_end()));
        }
        else
            inverted.push_back(regions[i]);
    }
}



/*************************************************
 * This function takes the names of three files (a reads file, a
 * chromosome file, and a dead zones file [possibly empty]), along
 * with a size for "deserts". The function returns (through
 * parameters) the "good" regions, the reads in the good region, and
 * the dead zones in each good region.
 */
void
LoadReadsByRegion(const char *regions_file, const char *reads_file, 
				  const char *deads_file, const size_t desert_size, 
				  const int VERBOSE,
				  vector<SimpleGenomicRegion> &regions, 
				  vector<vector<SimpleGenomicRegion> > &reads,
				  vector<vector<SimpleGenomicRegion> > &deads) 
{
  
    // get the regions
    if (VERBOSE)
        cerr << "[LOADING_DATA] chromosomes" << endl;
    vector<SimpleGenomicRegion> raw_regions;
    ReadBEDFile(regions_file, raw_regions);
    if (!check_sorted(raw_regions, true))
        SortGenomicRegion::sort_regions_collapse(raw_regions);
  
    // get the dead zones
    vector<SimpleGenomicRegion> raw_deads;
    if (deads_file) 
    {
        if (VERBOSE)
            cerr << "[LOADING_DATA] dead zones" << endl;
        ReadBEDFile(deads_file, raw_deads);
        if (!check_sorted(raw_deads, true))
            SortGenomicRegion::sort_regions_collapse(raw_deads);
    }
  
    // Get the deserts
    vector<SimpleGenomicRegion> deserts;
    get_deserts(raw_deads, desert_size, deserts);
  
    // Invert the deserts with respect to the regions
    invert_deserts(raw_regions, deserts, regions);
    if (!check_sorted(regions, true))
        SortGenomicRegion::sort_regions_collapse(regions);

    // Load in the file of reads
    if (VERBOSE)
        cerr << "[LOADING_DATA] reads" << endl;
    vector<SimpleGenomicRegion> raw_reads;
    ReadBEDFile(reads_file, raw_reads);
    if (!check_sorted(raw_reads, false))
        SortGenomicRegion::sort_regions(raw_reads);

    if (VERBOSE)
        cerr << "[LOADING_DATA] separating deserts" << endl;
    separate_regions(regions, raw_reads, reads);
    separate_regions(regions, raw_deads, deads);
}



void
LoadReadsByRegion(const char *regions_file, 
				  const char *reads_file_a, const char *reads_file_b, 
				  const char *deads_file, const size_t desert_size, 
				  const int VERBOSE,
				  vector<SimpleGenomicRegion> &regions,
				  vector<vector<SimpleGenomicRegion> > &reads_a,
				  vector<vector<SimpleGenomicRegion> > &reads_b,
				  vector<vector<SimpleGenomicRegion> > &deads) 
{
  
    // get the regions
    if (VERBOSE)
        cerr << "[LOADING_DATA] chromosomes" << endl;
    vector<SimpleGenomicRegion> raw_regions;
    ReadBEDFile(regions_file, raw_regions);
    if (!check_sorted(raw_regions, true))
        SortGenomicRegion::sort_regions_collapse(raw_regions);

  
    // get the dead zones
    vector<SimpleGenomicRegion> raw_deads;
    if (deads_file) 
    {
        if (VERBOSE)
            cerr << "[LOADING_DATA] dead zones" << endl;
        ReadBEDFile(deads_file, raw_deads);
        if (!check_sorted(raw_deads, true))
            SortGenomicRegion::sort_regions_collapse(raw_deads);
    }
  
    // Get the deserts
    vector<SimpleGenomicRegion> deserts;
    get_deserts(raw_deads, desert_size, deserts);
  
    // Invert the deserts with respect to the regions
    invert_deserts(raw_regions, deserts, regions);
    if (!check_sorted(regions, true))
        SortGenomicRegion::sort_regions_collapse(regions);

    // Load in the file of reads
    if (VERBOSE)
        cerr << "[LOADING_DATA] " << reads_file_a << endl;
    vector<SimpleGenomicRegion> raw_reads_a;
    ReadBEDFile(reads_file_a, raw_reads_a);
    if (!check_sorted(raw_reads_a, false))
        SortGenomicRegion::sort_regions(raw_reads_a);


    if (VERBOSE)
        cerr << "[LOADING_DATA] " << reads_file_b << endl;
    vector<SimpleGenomicRegion> raw_reads_b;
    ReadBEDFile(reads_file_b, raw_reads_b);
    if (!check_sorted(raw_reads_b, false))
        SortGenomicRegion::sort_regions(raw_reads_b);

    if (VERBOSE)
        cerr << "[LOADING_DATA] separating deserts" << endl;
    separate_regions(regions, raw_reads_a, reads_a);
    separate_regions(regions, raw_reads_b, reads_b);
    separate_regions(regions, raw_deads, deads);
}
