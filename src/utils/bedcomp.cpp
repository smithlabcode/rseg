#include <fstream>
#include <map>
#include <iomanip>
#include <cmath>
#include <string>

#include <popt.h>

#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::string;

int
main(int argc, const char **argv) {

  try {
    /* FILES */
//     const char *bedfilea = 0;
//     const char *bedfileb = 0;
//     const char *outfilename = 0;
//     const char *chromfile = 0;
		  string bedfilea, bedfileb, outfilename,chromfile;
		  

    double epsilon = 1e-40;
    int VERBOSE = 0;
  
    ////////////////////// COMMAND LINE OPTIONS /////////////////////////
//     static struct poptOption optionsTable[] = {
//       { "chrom", 'c', POPT_ARG_STRING, &chromfile, 0, 
// 	"chrom file" },
//       { "outfile", 'o', POPT_ARG_STRING, &outfilename, 0, 
// 	"output file (default: stdout)" },
//       { "epsilon", 'e', POPT_ARG_DOUBLE, &epsilon, 0, 
// 	"smallest probability value used in calculations" },
//       { "verbose", 'v', POPT_ARG_NONE, &VERBOSE, 0, 
// 	"print more run info" },
//       POPT_AUTOHELP POPT_TABLEEND
//     };

    OptionParser opt_parse("bedcomp", "A program for comparing BED format files",
			   "<File_A> <File_B>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfilename);
	opt_parse.add_opt("chrom", 'c', "chrom file", true, chromfile);
	opt_parse.add_opt("episeg", 'e', "smallest probability value used in calculations",
					  false, epsilon);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
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

  
//     /****************** GET COMMAND LINE ARGUMENTS ***************************/
//     poptContext optCon = poptGetContext("bedcomp", argc, argv, optionsTable, 0);
//     poptSetOtherOptionHelp(optCon, "<file_A> <file_B>");
//     if (argc < 2) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_SUCCESS;
//     }
//     char c;
//     if ((c = poptGetNextOpt(optCon)) < -1) {
//       cerr << "bedcomp: bad argument "
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
//       cerr << "bedcomp: leftover argument " << poptGetArg(optCon) << endl;
//       return EXIT_FAILURE;
//     }
//     poptFreeContext(optCon);
    /**********************************************************************/

    vector<SimpleGenomicRegion> chroms;
    ReadBEDFile(chromfile.c_str(), chroms);
    if (!check_sorted(chroms)) {
      cerr << "ERROR: " << chromfile << " not sorted" << endl;
      return EXIT_FAILURE;
    }
    double genome_size = 0;
    for (size_t i = 0; i < chroms.size(); ++i) {
      genome_size += chroms[i].get_width();
    }
    
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
    
    size_t intersection = 0;
    for (size_t i = 0; i < regions_c.size(); ++i)
      intersection += regions_c[i].get_width();

    size_t a_total = 0;
    for (size_t i = 0; i < regions_a.size(); ++i)
      a_total += regions_a[i].get_width();

    size_t b_total = 0;
    for (size_t i = 0; i < regions_b.size(); ++i)
      b_total += regions_b[i].get_width();
    
    size_t uunion = a_total + b_total - intersection; 

    if (VERBOSE)
      cerr << "genome size=" << genome_size << endl
	   << "a total    =" << a_total << endl
	   << "b total    =" << b_total << endl
	   << "intersct   =" << intersection << endl << endl;
    
    double p_x1 = a_total/genome_size;
    double p_y1 = 1 - p_x1;
    double p_x2 = b_total/genome_size;
    double p_y2 = 1 - p_x2;
    
    double p_xx = intersection/genome_size;
    double p_yy = 1 - uunion/genome_size;
    double p_xy = (uunion - b_total)/genome_size;
    double p_yx = (uunion - a_total)/genome_size;

    if (p_x1 < epsilon) p_x1 = epsilon;
    if (p_x2 < epsilon) p_x2 = epsilon;
    if (p_y1 < epsilon) p_y1 = epsilon;
    if (p_y2 < epsilon) p_y2 = epsilon;

    if (p_xx < epsilon) p_xx = epsilon;
    if (p_yy < epsilon) p_yy = epsilon;
    if (p_xy < epsilon) p_xy = epsilon;
    if (p_yx < epsilon) p_yx = epsilon;

    const double H_1 = -p_x1*log(p_x1) - p_y1*log(p_y1);
    const double H_2 = -p_x2*log(p_x2) - p_y2*log(p_y2);

    const double I = ((p_xx*(log(p_xx) - log(p_x1) - log(p_x2))) +
		      (p_xy*(log(p_xy) - log(p_x1) - log(p_y2))) +
		      (p_yy*(log(p_yy) - log(p_y1) - log(p_y2))) +
		      (p_yx*(log(p_yx) - log(p_y1) - log(p_x2))));
  
    if (VERBOSE) {
      cerr << "x_1 = " << bedfilea << endl
	   << "y_1 = " << "1 - " << bedfilea << endl
	   << "x_2 = " << bedfileb << endl
	   << "y_2 = " << "1 - " << bedfileb << endl << endl;
      
      cerr << "p_x1=" << p_x1 << endl
	   << "p_y1=" << p_y2 << endl
	   << "p_x2=" << p_x2 << endl
	   << "p_y2=" << p_y2 << endl << endl;
      
      cerr << "p_xx=" << p_xx << endl
	   << "p_xy=" << p_xy << endl
	   << "p_yx=" << p_yx << endl
	   << "p_yy=" << p_yy << endl << endl;
      
      cerr << "H(1)=" << H_1 << endl
	   << "H(2)=" << H_2 << endl
	   << "I(1,2)=" << I << endl << endl;
    }
    
    std::ostream* out = (!outfilename.empty()) ? new std::ofstream(outfilename.c_str()) : &std::cout;
    *out << "d_seg=" << 1 - 2*I/(H_1 + H_2) << endl
	 << "IntersectProp=" 
	 << static_cast<double>(a_total + b_total - 2*intersection)/
      (a_total + b_total) << endl;
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
