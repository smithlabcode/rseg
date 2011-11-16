#include <fstream>
#include <popt.h>

#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <GenomicProfile.hpp>

#include "OptionParser.hpp"

using std::ofstream;
using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::make_pair;
using std::pair;

template <class T> static void
BinReads(const vector<T> &reads,
	 const size_t region_start, const size_t region_end,
	 const size_t bin_size, vector<double> &bins) {
  bins.clear();
  size_t read_idx = 0;
  for (size_t i = region_start; i < region_end; i += bin_size) {
    size_t counts = 0;
    while (read_idx < reads.size() && reads[read_idx].get_start() < i + bin_size) {
      if (reads[read_idx].get_start() >= i)
	++counts;
      ++read_idx;
    }
    bins.push_back(counts);
  }
}


template <class T> static void
cluster_by_overlap(const vector<T> &mapped_locations, 
		   vector<vector<size_t> > &clusters, vector<T> &regions) {
  clusters = vector<vector<size_t> >(1, vector<size_t>(1, 0));
  regions.push_back(T(mapped_locations.front()));
  for (size_t i = 1; i < mapped_locations.size(); ++i) {
    if (!mapped_locations[i - 1].overlaps(mapped_locations[i])) {
      regions.back().set_end(mapped_locations[i - 1].get_end());
      regions.push_back(T(mapped_locations[i]));
      clusters.push_back(vector<size_t>());
    }
    clusters.back().push_back(i);
  }
  regions.back().set_end(mapped_locations.back().get_end());
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
		  string reads_file, outfile;
    int bin_size = 100;

//     ////////////////////// COMMAND LINE OPTIONS /////////////////////////
//     static struct poptOption optionsTable[] = {
//       { "output", 'o', POPT_ARG_STRING, &outfile, 0, 
// 	"output file (default: stdout)" },
//       { "bin", 'b', POPT_ARG_INT, &bin_size, 0, 
// 	"size of bins" },
//       POPT_AUTOHELP POPT_TABLEEND
//     };

    OptionParser opt_parse("readswig", "A program", "reads_file");
	opt_parse.add_opt("output", 'o', "output file name", false, outfile);
	opt_parse.add_opt("bin", 'b', "bin_size", false, bin_size);
	
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
	
	reads_file = leftover_args.front();


//     /****************** GET COMMAND LINE ARGUMENTS ***************************/
//     poptContext optCon = poptGetContext("readswig", argc, argv, optionsTable, 0);
//     poptSetOtherOptionHelp(optCon, "<reads>");
//     if (argc < 2) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_SUCCESS;
//     }
//     char c;
//     if ((c = poptGetNextOpt(optCon)) < -1) {
//       cerr << "readswig: bad argument "
// 	   << poptBadOption(optCon, POPT_BADOPTION_NOALIAS) << ": "
// 	   << poptStrerror(c) << endl;
//       return EXIT_FAILURE;
//     }
//     if (poptPeekArg(optCon) == 0) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_FAILURE;
//     }
//     else reads_file = poptGetArg(optCon);
//     if (poptPeekArg(optCon) != 0) {
//       cerr << "readswig: leftover argument " << poptGetArg(optCon) << endl;
//       return EXIT_FAILURE;
//     }
//     poptFreeContext(optCon);
//     /**********************************************************************/

    vector<SimpleGenomicRegion> reads;
    ReadBEDFile(reads_file.c_str(), reads);
    assert(check_sorted(reads));
    
    vector<SimpleGenomicRegion> regions;
    vector<vector<size_t> > clusters;
    cluster_by_overlap(reads, clusters, regions);
    
    ostream* out = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
    for (size_t i = 0; i < regions.size(); ++i) {
      vector<double> values;
      vector<SimpleGenomicRegion> clust;
      for (size_t j = 0; j < clusters[i].size(); ++j)
	clust.push_back(reads[clusters[i][j]]);
      BinReads(clust, regions[i].get_start(), regions[i].get_end(), bin_size, values);
      GenomicProfile prof(regions[i], values, bin_size);
      *out << prof << endl;
    }
    if (out != &cout) delete out;
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
