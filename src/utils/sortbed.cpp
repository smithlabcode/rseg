/*    sortbed: a program for sorting BED format files
 *    Copyright (C) 2008 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

typedef GenomicRegion* GenomicRegionPointer;

struct region_pointer_less {
  bool operator()(const GenomicRegionPointer a, 
		  const GenomicRegionPointer b) const {
    return (*a) < (*b);
  }
};

void
sort_regions(vector<GenomicRegion> &regions) {
  vector<GenomicRegionPointer> sorter;
  for (vector<GenomicRegion>::iterator i = regions.begin(); 
       i != regions.end(); ++i) sorter.push_back(&(*i));
  sort(sorter.begin(), sorter.end(), region_pointer_less());
  
  vector<GenomicRegion> r;
  r.reserve(regions.size());
  for (vector<GenomicRegionPointer>::const_iterator i(sorter.begin());
       i != sorter.end(); ++i)
    r.push_back(*(*i));
  r.swap(regions);
}

void
sort_regions_collapse_chrom(vector<GenomicRegion> &regions) {
  static const string FAKE_NAME("X");
  const string chrom(regions.front().get_chrom());

  vector<pair<size_t, bool> > boundaries;
  for (size_t i = 0; i < regions.size(); ++i) {
    boundaries.push_back(make_pair(regions[i].get_start(), false));
    boundaries.push_back(make_pair(regions[i].get_end(), true));
  }
  regions.clear();
  sort(boundaries.begin(), boundaries.end());

  GenomicRegion holder(chrom, 0, 0, FAKE_NAME, 0, '+');
  size_t count = 0;
  for (size_t i = 0; i < boundaries.size(); ++i)
    if (boundaries[i].second) {
      --count;
      if (count == 0) {
	holder.set_end(boundaries[i].first);
	regions.push_back(holder);
      }
    }
    else {
      if (count == 0) 
	holder.set_start(boundaries[i].first);
      ++count;
    }
}

void
sort_regions_collapse(vector<GenomicRegion> &regions) {
  vector<vector<GenomicRegion> > separated_by_chrom;
  separate_chromosomes(regions, separated_by_chrom);
  regions.clear();
  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    sort_regions_collapse_chrom(separated_by_chrom[i]);
    regions.insert(regions.end(),
		   separated_by_chrom[i].begin(), 
		   separated_by_chrom[i].end());
    separated_by_chrom[i].clear();
  }
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool collapse_regions = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("sortbed", "A program for sorting BED format files",
			   "<bed-format-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("collapse", 'c', "Collapse the BED intervals", 
		      false , collapse_regions);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file_name = leftover_args.front();
    /**********************************************************************/
    
    vector<GenomicRegion> regions;
    ReadBEDFile(input_file_name, regions);

    if (collapse_regions)
      sort_regions_collapse(regions);
    else if (!check_sorted(regions))
      sort_regions(regions);
    
    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    copy(regions.begin(), regions.end(), 
	 std::ostream_iterator<GenomicRegion>(*out, "\n"));
    if (out != &std::cout) delete out;
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
