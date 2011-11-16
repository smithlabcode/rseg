
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

#include "string"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::string;

int
main(int argc, const char **argv) {
  
  try {
  
    /* FILES */
//     const char *bedfilea = 0;
//     const char *bedfileb = 0;

    ////////////////////// COMMAND LINE OPTIONS /////////////////////////
//     static struct poptOption optionsTable[] = {
//       POPT_AUTOHELP POPT_TABLEEND
//     };
    
//     /****************** GET COMMAND LINE ARGUMENTS ***************************/
//     poptContext optCon = poptGetContext("bedintersect", argc, argv, optionsTable, 0);
//     poptSetOtherOptionHelp(optCon, "<regions_a> <regions_b>");
//     if (argc < 2) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_SUCCESS;
//     }
//     char c;
//     if ((c = poptGetNextOpt(optCon)) < -1) {
//       cerr << "bedintersect: bad argument "
// 	   << poptBadOption(optCon, POPT_BADOPTION_NOALIAS) << ": "
// 	   << poptStrerror(c) << endl;
//       return EXIT_FAILURE;
//     }
//     if (poptPeekArg(optCon) == 0) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_FAILURE;
//     }
//     else bedfilea = poptGetArg(optCon);
//     if (poptPeekArg(optCon) == 0) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_FAILURE;
//     }
//     else bedfileb = poptGetArg(optCon);
//     if (poptPeekArg(optCon) != 0) {
//       cerr << "bedintersect: leftover argument " << poptGetArg(optCon) << endl;
//       return EXIT_FAILURE;
//     }
//     poptFreeContext(optCon);

    OptionParser opt_parse("bedintersect", "A program for intersecting two bed files",
			   "<File_A> <File_B>");
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
	if (leftover_args.size() < 2)
	{
			cerr << "Need two file names" << endl;
			return EXIT_SUCCESS;
	}
	
    const string bedfilea = leftover_args[0];
	const string bedfileb = leftover_args[1];


    /**********************************************************************/

    vector<SimpleGenomicRegion> regions_a;
    ReadBEDFile(bedfilea.c_str(), regions_a);
    if (!check_sorted(regions_a)) {
      cerr << "ERROR: " << bedfilea << " not sorted" << endl;
      return EXIT_FAILURE;
    }
    vector<SimpleGenomicRegion> regions_b;
    ReadBEDFile(bedfileb.c_str(), regions_b);
    if (!check_sorted(regions_b)) {
      cerr << "ERROR: " << bedfileb << " not sorted" << endl;
      return EXIT_FAILURE;
    }
    vector<SimpleGenomicRegion> regions_c;
    genomic_region_intersection_by_base(regions_a, regions_b, regions_c);
    
    copy(regions_c.begin(), regions_c.end(), 
	 std::ostream_iterator<SimpleGenomicRegion>(cout, "\n"));
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
