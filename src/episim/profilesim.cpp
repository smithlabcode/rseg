#include <fstream>
#include <map>
#include <iomanip>
#include <numeric>

#include <popt.h>

#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "GenomicProfile.hpp"
#include "Distro.hpp"

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
using std::setprecision;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*
static void
simulate_profile(gsl_rng *rng, const SimpleGenomicRegion &chrom,
		 const vector<SimpleGenomicRegion> &domains,
		 const size_t step, 
		 const double fg_gamma_k, const double fg_gamma_theta,
		 const double bg_gamma_k, const double bg_gamma_theta,
		 vector<double> &profile) {
  
  const size_t chrom_size = chrom.get_width();

  size_t start_pos = domains.front().get_start();
  size_t end_pos = domains.front().get_end();

  size_t i = 0, j = 0;
  while (i < chrom_size && j < domains.size()) {
    if (i > start_pos) {
      while (i < end_pos) {
	profile.push_back(gsl_ran_gamma(rng, fg_gamma_k, fg_gamma_theta));
	i += step;
      }
      ++j;
      if (j < domains.size()) {
	start_pos = domains[j].get_start();
	end_pos = domains[j].get_end();
      }
    }
    else while (i <= start_pos) {
	profile.push_back(gsl_ran_gamma(rng, bg_gamma_k, bg_gamma_theta));
	i += step;
      }
  }
  while (i < chrom_size) {
    profile.push_back(gsl_ran_gamma(rng, bg_gamma_k, bg_gamma_theta));
    i += step;
  }
}


static void
simulate_profile(gsl_rng *rng, 
		 const SimpleGenomicRegion &chrom, const size_t step, 
		 const double bg_gamma_k, const double bg_gamma_theta,
		 vector<double> &profile) {
  const size_t chrom_size = chrom.get_width();
  for (size_t i = 0; i < chrom_size; i += step)
    profile.push_back(gsl_ran_gamma(rng, bg_gamma_k, bg_gamma_theta));
}



static void
simulate_profile(gsl_rng *rng, 
		 const vector<SimpleGenomicRegion> &chroms, 
		 const vector<SimpleGenomicRegion> &domains, 
		 const size_t step, 
		 const double fg_gamma_k, const double fg_gamma_theta,
		 const double bg_gamma_k, const double bg_gamma_theta,
		 vector<GenomicProfile> &profile) {
  
  vector<vector<SimpleGenomicRegion> > sep_domains;
  separate_chromosomes(domains, sep_domains);
  
  map<string, size_t> chrom_dom;
  for (size_t i = 0; i < sep_domains.size(); ++i)
    chrom_dom[sep_domains[i].front().get_chrom()] = i;
  
  vector<vector<double> > values(chroms.size());
  for (size_t i = 0; i < chroms.size(); ++i) {
    const map<string, size_t>::const_iterator d = 
      chrom_dom.find(chroms[i].get_chrom());
    if (d != chrom_dom.end())
      simulate_profile(rng, chroms[i], sep_domains[d->second], step, 
		       fg_gamma_k, fg_gamma_theta,
		       bg_gamma_k, bg_gamma_theta, values[i]);
    else simulate_profile(rng, chroms[i], step, 
			  bg_gamma_k, bg_gamma_theta, values[i]);
  }
  
  double total = 0;
  for (size_t i = 0; i < values.size(); ++i)
    total += std::accumulate(values[i].begin(), values[i].end(), 0.0);
  
  for (size_t i = 0; i < values.size(); ++i) {
    std::transform(values[i].begin(), values[i].end(), values[i].begin(),
		   std::bind2nd(std::divides<double>(), total));
    profile.push_back(GenomicProfile(chroms[i], values[i], step));
    values.clear();
  }
}



int
main(int argc, const char **argv) {

  const char *program_description =
    "The purpose of this program is to take a set of genomic domains and generate\n"
    "a profile with a particular mean level inside the domains, and a different\n"
    "mean level outside the domains.";

  const char *chroms_file = 0;
  const char *domains_file = 0;
  const char *outfile = 0;

  double fg_gamma_k = 2.0;
  double fg_gamma_theta = 1.0;
  
  double bg_gamma_k = 1.0;
  double bg_gamma_theta = 1.0;
  
  int step = 5;
  
  int VERBOSE = 0;
  size_t seed = std::numeric_limits<size_t>::max();

  const char *profile_name_input = 0;
  
  ////////////////////// COMMAND LINE OPTIONS /////////////////////////
  static struct poptOption optionsTable[] = {
    { "domains", 'd', POPT_ARG_STRING, &domains_file, 0, 
      "file of domains" },
    { "fg-k", 'K', POPT_ARG_DOUBLE, &fg_gamma_k, 0, 
      "foreground gamma a parameter" },
    { "bg-k", 'k', POPT_ARG_DOUBLE, &bg_gamma_k, 0, 
      "background gamma a parameter" },
    { "fg-theta", 'T', POPT_ARG_DOUBLE, &fg_gamma_theta, 0, 
      "foreground gamma b parameter" },
    { "bg-theta", 't', POPT_ARG_DOUBLE, &bg_gamma_theta, 0, 
      "background gamma b parameter" },
    { "step", 's', POPT_ARG_INT, &step, 0, 
      "step size" },
    { "out", 'o', POPT_ARG_STRING, &outfile, 0, 
      "output file name" },
    { "verbose", 'v', POPT_ARG_NONE, &VERBOSE, 0, 
      "print more run information" },
    { "seed", 'S', POPT_ARG_INT, &seed, 0, 
      "random number seed" },
    { "name", 'N', POPT_ARG_STRING, &profile_name_input, 0, 
      "name of profile" },
    POPT_AUTOHELP POPT_TABLEEND
  };
  
  // ****************** GET COMMAND LINE ARGUMENTS ***************************
  poptContext optCon = poptGetContext("profilesim", argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "<chroms-file>");
  if (argc < 2) {
    poptPrintHelp(optCon, stderr, 0);
    cerr << endl << program_description << endl << endl;
    return EXIT_SUCCESS;
  }
  char c;
  if ((c = poptGetNextOpt(optCon)) < -1) {
    cerr << "profilesim: bad argument "
	 << poptBadOption(optCon, POPT_BADOPTION_NOALIAS) << ": "
	 << poptStrerror(c) << endl;
    return EXIT_FAILURE;
  }
  if (poptPeekArg(optCon) == 0) {
    poptPrintHelp(optCon, stderr, 0);
    return EXIT_FAILURE;
  }
  else chroms_file = poptGetArg(optCon);
  if (poptPeekArg(optCon) != 0) {
    cerr << "profilesim: leftover argument " 
	 << poptGetArg(optCon) << endl;
    return EXIT_FAILURE;
  }
  poptFreeContext(optCon);
  // **********************************************************************
  
  try {

    if (VERBOSE)
      cerr << "fg_mean/bg_mean=" 
	   << fg_gamma_k*fg_gamma_theta/
	(bg_gamma_k*bg_gamma_theta) << endl;

    // Setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    if (seed != std::numeric_limits<size_t>::max())
      srand(seed);
    else {
      const int random_number_seed = time(0) + getpid();
      if (VERBOSE)
	cerr << "random number seed=" << random_number_seed << endl;
      srand(random_number_seed);
    }
    gsl_rng_set(rng, rand());
    
    if (VERBOSE)
      cerr << "reading chroms" << endl;
    vector<SimpleGenomicRegion> chroms;
    ReadBEDFile(chroms_file, chroms);
    check_sorted(chroms, true);

    if (VERBOSE)
      cerr << "reading domains" << endl;
    vector<SimpleGenomicRegion> domains;
    ReadBEDFile(domains_file, domains);
    check_sorted(domains, true);
    
    if (VERBOSE)
      cerr << "simulating" << endl;
    vector<GenomicProfile> profile;
    simulate_profile(rng, chroms, domains, step, 
		     fg_gamma_k, fg_gamma_theta, 
		     bg_gamma_k, bg_gamma_theta, profile);
    const string profile_name((profile_name_input) ? profile_name_input : "");
    WriteWIGFile(outfile, profile_name, profile);
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
