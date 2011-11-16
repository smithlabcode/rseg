#include <fstream>
#include <map>
#include <iomanip>
#include <numeric>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "smithlab_utils.hpp"
#include "GenomicProfile.hpp"
#include "OptionParser.hpp"

using std::ofstream;
using std::string;
using std::vector;
using std::map;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::pair;
using std::setw;
using std::make_pair;
using std::max;
using std::min;
using std::ios_base;
using std::accumulate;
using std::ptr_fun;

/*
void
simulate_from_profile(gsl_rng *rng, const char *profile_file, 
		      const size_t n_reads,
		      vector<SimpleGenomicRegion> &reads) {
  
  vector<GenomicProfile> profiles;
  ReadWIGFile(profile_file, profiles);
  
  for (size_t i = 0; i < profiles.size(); ++i) {
    const string chrom_name(profiles[i].get_region().get_chrom());
    const size_t chrom_start = profiles[i].get_region().get_start();
    const size_t step = profiles[i].get_step();
    const size_t profile_size = profiles[i].get_size();
    
    vector<unsigned int> read_counts(profile_size);
    gsl_ran_multinomial(rng, profile_size, n_reads, 
			&profiles[i][0], &read_counts.front());
    for (size_t j = 0; j < profile_size; ++j) {
      for (size_t k = 0; k < read_counts[j]; ++k) {
	const size_t start = (chrom_start + j*step + 
			      gsl_rng_uniform_int(rng, step));
	reads.push_back(SimpleGenomicRegion(chrom_name, start, start + 1));
      }
    }
  }
}

int
main(int argc, const char **argv) {
  
  try {

    const char *profile_file = 0;
  
    const char *outfile = 0;
    size_t seed = std::numeric_limits<size_t>::max();

    size_t n_reads = 0;
    int VERBOSE = 0;
  
    ////////////////////// COMMAND LINE OPTIONS /////////////////////////
    static struct poptOption optionsTable[] = {
      { "reads", 'r', POPT_ARG_INT, &n_reads, 0, 
	"number of reads to sequence" },
      { "out", 'o', POPT_ARG_STRING, &outfile, 0, 
	"output file name" },
      { "verbose", 'v', POPT_ARG_NONE, &VERBOSE, 0, 
	"print more run information" },
      { "seed", 'S', POPT_ARG_INT, &seed, 0, 
	"rng seed" },
      POPT_AUTOHELP POPT_TABLEEND
    };
  
    // ****************** GET COMMAND LINE ARGUMENTS ***************************
    poptContext optCon = poptGetContext("histmodsim", argc, argv, optionsTable, 0);
    poptSetOtherOptionHelp(optCon, "<profile.wig>");
    if (argc < 2) {
      poptPrintHelp(optCon, stderr, 0);
      return EXIT_SUCCESS;
    }
    char c;
    if ((c = poptGetNextOpt(optCon)) < -1) {
      cerr << "histmodsim: bad argument "
	   << poptBadOption(optCon, POPT_BADOPTION_NOALIAS) << ": "
	   << poptStrerror(c) << endl;
      return EXIT_FAILURE;
    }
    if (poptPeekArg(optCon) == 0) {
      poptPrintHelp(optCon, stderr, 0);
      return EXIT_FAILURE;
    }
    else profile_file = poptGetArg(optCon);
    if (poptPeekArg(optCon) != 0) {
      cerr << "histmodsim: leftover argument " 
	   << poptGetArg(optCon) << endl;
      return EXIT_FAILURE;
    }
    poptFreeContext(optCon);
    // **********************************************************************
  
    // Setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    if (seed == std::numeric_limits<size_t>::max())
      seed = time(0) + getpid();
    if (VERBOSE)
      cerr << "random number seed=" << seed << endl;
    srand(seed);
    gsl_rng_set(rng, rand());
    
    // get the chromosomes
    vector<SimpleGenomicRegion> reads;
    simulate_from_profile(rng, profile_file, n_reads, reads);
    sort(reads.begin(), reads.end());
    
    if (VERBOSE) 
      cerr << "writing reads" << endl;
    std::ostream* out = (outfile) ? new std::ofstream(outfile) : &cout;
    copy(reads.begin(), reads.end(), 
	 std::ostream_iterator<SimpleGenomicRegion>(*out, "\n"));
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
*/

int 
main()
{
}
