#include <fstream>
#include <map>
#include <iomanip>
#include <numeric>

// #include <popt.h>

#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <ReadCounts.hpp>
#include <LoadReadsByRegion.hpp>
#include <OptionParser.hpp>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_histogram.h>

using std::ofstream;
using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
/*
void
reads_hist(vector<vector<double> > &counts,
		   int min_val, int max_val, gsl_histogram *h) {
		gsl_histogram_set_ranges_uniform(h, min_val, max_val);
		for (size_t i = 0; i < counts.size(); ++i)
				for (size_t j = 0; j < counts[i].size(); ++j)
						gsl_histogram_increment(h, counts[i][j]);
}


int
main(int argc, const char **argv) {

		// file names
		string reads_file = 0;
		string deads_file = 0;
		string chroms_file = 0;
		string outfile = 0;
  
		// flags
		bool VERBOSE = false;
  
		size_t desert_size = std::numeric_limits<int>::max()/2;
		size_t bin_size = 1000;
  
		// commandline option parsing 
		OptionParser opt_parse("readcounthist",
							   "Estimate reads density",
							   "input-file-name");
		opt_parse.add_opt("chrom", 'c', "BED file with sizes of chromosomes",
						  true, chroms_file);
		opt_parse.add_opt("bin", 'b', "Bin size",
						  false, bin_size);
		opt_parse.add_opt("deadzone", 'd', "Filename of deadzones",
						  false, deads_file);
		opt_parse.add_opt("desert", 'D', "Desert size",
						  false, desert_size);
		opt_parse.add_opt("output", 'o', "Output file name",
						  false, outfile);
		opt_parse.add_opt("verbose", 'v', "Print more information",
						  false, VERBOSE);
		vector<string> leftover_args;
		try
		{
				opt_parse.parse(argc, argv, leftover_args);
		}
		catch(OptionParserException e)
		{
				cerr << e.what() << endl;
				return EXIT_FAILURE;
		}
		
		if (argc == 1 || opt_parse.help_requested())
		{
				cerr << opt_parse.help_message() << endl;
				return EXIT_SUCCESS;
		}
		if (opt_parse.about_requested())
		{
				cerr << opt_parse.about_message() << endl;
				return EXIT_SUCCESS;
		}
		if (opt_parse.option_missing())
		{
				cerr << opt_parse.option_missing_message() << endl;
				return EXIT_SUCCESS;
		}
		if (leftover_args.empty())
		{
				cerr << opt_parse.help_message() << endl;
				return EXIT_SUCCESS;
		}
		reads_file = leftover_args.front();
		


		try {
    
				vector<SimpleGenomicRegion> regions;
				vector<vector<SimpleGenomicRegion> > reads;
				vector<vector<SimpleGenomicRegion> > deads;
				LoadReadsByRegion(chroms_file.c_str(), reads_file.c_str(), deads_file.c_str(),
								  desert_size, VERBOSE,
								  regions, reads, deads);
    
    
				// Obtain the binned reads
				if (VERBOSE)
						cerr << "[preparing data] binning reads" << endl;
				vector<vector<double> > read_bins;
				vector<vector<SimpleGenomicRegion> > bin_boundaries;
				BinReadsExtendOverDeadZones(reads, deads, regions,
											bin_size, bin_boundaries, 
											read_bins);

				// make the histogram
				gsl_histogram *h = gsl_histogram_alloc(bin_size);
				reads_hist(read_bins, 0, bin_size, h);
    
				FILE *f = (!outfile.empty()) ? fopen(outfile.c_str(), "w") : stdout;
				gsl_histogram_fprintf(f, h, "%g", "%g");
				if (f != stdout) fclose(f);
				gsl_histogram_free(h);
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
*/

int 
main()
{
}
