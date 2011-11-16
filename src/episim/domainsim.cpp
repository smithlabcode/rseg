/*
 * Copyright (C) 2008 Cold Spring Harbor Laboratory
 *                    Andrew D Smith
 * Author: Andrew D. Smith
 *
 * This file is part of RMAP
 *
 * RMAP is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * RMAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * To receive a copy of the GNU General Public License, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301 USA
 */

#include <fstream>
#include <map>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <functional>
#include <string>


#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "Distro.hpp"
#include "RNG.hpp"
#include "OptionParser.hpp"


using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;
using std::max;
using std::multiplies;
using std::bind2nd;

/*
static void
split_first_interdomain_size(const size_t seed,
			     vector<double> &inter_domain_sizes) {
  assert(!inter_domain_sizes.empty());
  Runif runif(seed);
  const double last_size = runif.runif(0.0, inter_domain_sizes.front());
  if (last_size > 0.0) {
    inter_domain_sizes.front() -= last_size;
    inter_domain_sizes.push_back(last_size);
  }
}


static void
build_domains(const string chrom_name, 
	      const vector<double> &domain_sizes, 
	      const vector<double> &inter_domain_sizes, 
	      vector<SimpleGenomicRegion> &domains) {
  
  size_t curr_pos = 0;
  for (size_t i = 0; i < domain_sizes.size(); ++i) {
    const size_t start = curr_pos + static_cast<size_t>(inter_domain_sizes[i]);
    const size_t end = start + static_cast<size_t>(domain_sizes[i]);
    domains.push_back(SimpleGenomicRegion(chrom_name, start, end));
    curr_pos = end;
  }
}


static void
scale_sizes(const double scale_to, vector<double> &sizes) {
  const double unscaled_total = 
    accumulate(sizes.begin(), sizes.end(), 0.0);
  const double scale = scale_to/unscaled_total;
  vector<double> tmp_sizes;
  transform(sizes.begin(), sizes.end(),
	    back_inserter(tmp_sizes), bind2nd(multiplies<double>(), scale));
  sizes.swap(tmp_sizes);
}

int
main(int argc, const char **argv) {

  try {

    const char *program_description =
      "The purpose of this program is to randomly generate a set\n"
      "of genomic regions, all contained inside a given set of regions.";
  
//     const char *outfile = 0;
//     const char *chrom_name = "chr1";
    string outfile, chrom_name("chr1");

    size_t chrom_size = 1000000;
    size_t n_domains = 0;
    double coverage = 0.5;
    size_t min_domain_size = 100;
  
//     const char *size_distro = "exp,1";
//     const char *inter_distro = "exp,1";
    string size_distro("exp,1"), inter_distro("exp,1");

    bool VERBOSE = 0;
    size_t seed = std::numeric_limits<size_t>::max();
    
    ////////////////////// COMMAND LINE OPTIONS /////////////////////////
//     static struct poptOption optionsTable[] = {
//       { "domains", 'n', POPT_ARG_INT, &n_domains, 0, 
// 	"number of domains to simulate" },
//       { "min", 'm', POPT_ARG_INT, &min_domain_size, 0, 
// 	"minimum size of domains" },
//       { "chrom", 'C', POPT_ARG_INT, &chrom_size, 0, 
// 	"size of simulated chrom" },
//       { "coverage", 'c', POPT_ARG_DOUBLE, &coverage, 0, 
// 	"percent of genome covered" },
//       { "sizes", 's', POPT_ARG_STRING, &size_distro, 0, 
// 	"distribution of domain sizes (with params)" },
//       { "out", 'o', POPT_ARG_STRING, &outfile, 0, 
// 	"output file name" },
//       { "verbose", 'v', POPT_ARG_NONE, &VERBOSE, 0, 
// 	"print more run information" },
//       { "seed", 'S', POPT_ARG_INT, &seed, 0, 
// 	"random number seed" },
//       POPT_AUTOHELP POPT_TABLEEND
//     };

    OptionParser opt_parse("domainsim",
                           " The purpose of this program is to randomly generate a set\n" +
                           "of genomic regions, all contained inside a given set of regions.",
                           "chroms_file");
    opt_parse.add_opt("domains", 'n', "Number of domains to simulate", false, n_domains);
    opt_parse.add_opt("min", 'm', "Minimal size of simulated domains", false, min_domain_size);
    opt_parse.add_opt("chrom", 'C', "Size of simulated chromosomes", false, chrom_size);
    opt_parse.add_opt("coverage", 'c', "Percent of genome covered", false, coverage);
    opt_parse.add_opt("sizes", 's', "distribution of domain sizes (with params)", false, size_distro);
    opt_parse.add_opt("out", 'o', "output file name", false, outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", false, VERBOSE);
    opt_parse.add_opt("seed", 'S', "random number seed", false, seed);

//     ****************** GET COMMAND LINE ARGUMENTS ***************************
//     poptContext optCon = poptGetContext("domainsim", argc, argv, optionsTable, 0);
//     poptSetOtherOptionHelp(optCon, "<chroms-file>");
//     if (argc < 2) {
//       poptPrintHelp(optCon, stderr, 0);
//       cerr << endl << program_description << endl << endl;
//       return EXIT_SUCCESS;
//     }
//     char c;
//     if ((c = poptGetNextOpt(optCon)) < -1) {
//       cerr << "domainsim: bad argument "
// 	   << poptBadOption(optCon, POPT_BADOPTION_NOALIAS) << ": "
// 	   << poptStrerror(c) << endl;
//       return EXIT_FAILURE;
//     }
//     if (poptPeekArg(optCon) != 0) {
//       cerr << "domainsim: leftover argument " 
// 	   << poptGetArg(optCon) << endl;
//       return EXIT_FAILURE;
//     }
//     poptFreeContext(optCon);

    
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


    if (!Distro::has_params(size_distro)) {
      cerr << "domainsim: need parameters for \"" << size_distro << "\"" << endl;
      return EXIT_FAILURE;
    }
    if (!Distro::has_params(inter_distro)) {
      cerr << "domainsim: need parameters for \"" << inter_distro << "\"" << endl;
      return EXIT_FAILURE;
    }


    
    if (VERBOSE)
      cerr << "sampling domain sizes" << endl;
    
    if (seed == std::numeric_limits<size_t>::max())
      seed = time(0) + getpid();
    
    Distro d(size_distro);
    d.seed(seed);
    vector<double> domain_sizes;
    for (size_t i = 0; i < n_domains; ++i)
      domain_sizes.push_back(d.sample());
    scale_sizes(coverage*chrom_size, domain_sizes);
    // std::generate_n(back_inserter(unscaled_domain_sizes), n_domains, d);
    vector<double> tmp_dom_sizes;
    for (size_t i = 0; i < n_domains; ++i)
      if (domain_sizes[i] > min_domain_size)
	tmp_dom_sizes.push_back(domain_sizes[i]);
    domain_sizes.swap(tmp_dom_sizes);
    scale_sizes(coverage*chrom_size, domain_sizes);
    
    Distro inter_d(inter_distro);
    inter_d.seed(seed);
    vector<double> inter_sizes;
    for (size_t i = 0; i < n_domains; ++i)
      inter_sizes.push_back(inter_d.sample());
    scale_sizes((1.0 - coverage)*chrom_size, inter_sizes);
    vector<double> tmp_inter_sizes;
    for (size_t i = 0; i < n_domains; ++i)
      if (inter_sizes[i] > min_domain_size)
	tmp_inter_sizes.push_back(inter_sizes[i]);
    inter_sizes.swap(tmp_inter_sizes);
    scale_sizes((1.0 - coverage)*chrom_size, inter_sizes);
    
    split_first_interdomain_size(seed, inter_sizes);
    vector<SimpleGenomicRegion> domains;
    build_domains(chrom_name, domain_sizes, inter_sizes, domains);
    
    if (VERBOSE) {
      cerr << "chrom size             = " << chrom_size << endl
	   << "dom size               = " << chrom_size*coverage << endl
	   << "inter size             = " << chrom_size*(1 - coverage) << endl
	   << "n domains              = " << domain_sizes.size() << endl
	   << "mean domain size       = " << (chrom_size*coverage)/n_domains << endl
	   << "mean inter-domain dist = " << (chrom_size*(1 - coverage))/n_domains << endl;
    }
    cout << chrom_name << "\t" << 0 << "\t" << chrom_size << endl;
    
    std::ostream* out = (outfile) ? new std::ofstream(outfile) : &cout;
    *out << "track name=domainsim" << endl;
    copy(domains.begin(), domains.end(), 
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
