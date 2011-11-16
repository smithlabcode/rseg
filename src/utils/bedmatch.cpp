
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

#include <fstream>

using std::ofstream;
using std::string;
using std::vector;
using std::max;
using std::min;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::pair;
using std::make_pair;

template <class T> size_t
intersection_size(const T &a, const T &b) {
//   const size_t high = min(a.get_end(), b.get_end());
//   const size_t low = max(a.get_start(), b.get_start());
//   return ((high <= low) ? 0 : high - low);
  return 
    (max(a.get_start(), b.get_start()) - min(a.get_start(), b.get_start())) +
    (max(a.get_end(), b.get_end()) - min(a.get_end(), b.get_end()));
}

struct Dat {
  size_t id;
  double prev;
  size_t prev_id;
  double mid;
  size_t mid_id;
  double next,;
  size_t next_id;
  double skip;
  size_t skip_id_a, skip_id_b;
  bool type;
  Dat() : id(FAKE), 
	  prev(0), prev_id(NONE), 
	  mid(0), mid_id(NONE), 
	  next(0), next_id(NONE),
	  skip(0),
	  skip_id_a(NONE),
	  skip_id_b(NONE), 
	  type(false) {}
  string 
  tostring() const {
    std::ostringstream ss;
    if (id == FAKE)
      ss << "F";
    else 
      ss << id;
    ss << "\t";
    
    if (prev_id == NONE)
      ss << "N";
    else if (prev_id == FAKE)
      ss << "F";
    else
      ss << prev_id;
    ss << "," << prev << "\t";

    if (mid_id == NONE)
      ss << "N";
    else if (mid_id == FAKE)
      ss << "F";
    else
      ss << mid_id;
    ss << "," << mid << "\t";
    
    if (next_id == NONE)
      ss << "N";
    else if (next_id == FAKE)
      ss << "F";
    else
      ss << next_id;
    ss << "," << next << "\t";
    
    if (skip_id_a == NONE)
      ss << "N";
    else if (skip_id_a == FAKE)
      ss << "F";
    else
      ss << skip_id_a;
    ss << ",";
    if (skip_id_b == NONE)
      ss << "N";
    else if (skip_id_b == FAKE)
      ss << "F";
    else
      ss << skip_id_b;
    ss << "," << skip << "\t" << type;
    return ss.str();
  }

  static const size_t NONE = static_cast<size_t>(-1);
  static const size_t FAKE = static_cast<size_t>(-2);
};


ostream&
operator<<(ostream &ss, const Dat &d) {
  return ss << d.tostring();
}

void
get_triple_scores(const vector<pair<size_t, SimpleGenomicRegion> > &a,
		  const vector<pair<size_t, SimpleGenomicRegion> > &b,
		  vector<Dat> &scores) {
  
  size_t j = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    scores.push_back(Dat());
    scores.back().id = a[i].first;
    
    double overlap_prev = 0;
    if (j < b.size() && 
	b[j].second.get_start() < a[i].second.get_start() &&
	b[j].second.get_end() > a[i].second.get_start() && 
	b[j].second.get_end() <= a[i].second.get_end()) {
      overlap_prev = intersection_size(b[j].second, a[i].second);
      scores.back().prev = overlap_prev;
      scores.back().prev_id = b[j].first;
      ++j;
    }
    
    while (j < b.size() && a[i].second.contains(b[j].second)) {
      const double overlap = intersection_size(b[j].second, a[i].second);//b[j].second.get_width();
      if (overlap > scores.back().mid) {
	scores.back().mid = overlap;
	scores.back().mid_id = b[j].first;
      }
      ++j;
    }

    if (j < b.size() && b[j].second.contains(a[i].second)) {
      scores.back().mid = intersection_size(b[j].second, a[i].second);//a[i].second.get_width();
      scores.back().type = true;
      scores.back().mid_id = b[j].first;
    }

    double overlap_next = 0;
    if (j < b.size() && 
	a[i].second.get_start() < b[j].second.get_start() &&
	a[i].second.get_end() > b[j].second.get_start() &&
	a[i].second.get_end() <= b[j].second.get_end()) {
      overlap_next = intersection_size(b[j].second, a[i].second);//a[i].second.get_end() - b[j].second.get_start();
      scores.back().next = overlap_next;
      scores.back().next_id = b[j].first;
    }
  }
}


void
assign_parts(const vector<Dat> &scores, 
	     double final_skip,
	     size_t final_a,
	     size_t final_b,
	     vector<pair<size_t, size_t> > &path) {
  
  vector<vector<double> > v(scores.size(), 
			    vector<double>(3, -numeric_limits<double>::max()));
  vector<vector<size_t> > trace(scores.size(), vector<size_t>(4));

  v.front()[0] = scores.front().prev;
  v.front()[1] = scores.front().mid + scores.front().skip;
  v.front()[2] = scores.front().next + scores.front().skip;

  for (size_t i = 1; i < scores.size(); ++i) {
    /* handle case of matching with prev:
     * Only 2 cases: [\_,\_] OR [_|_,\_]
     */
    // CASE 1: [\_,\_]
    double prev_prev = scores[i].prev + v[i - 1][0];
    // CASE 2: [_|_,\_]
    double mid_prev = scores[i].prev + v[i - 1][1];

    if (prev_prev > mid_prev) {
      v[i][0] = prev_prev;
      trace[i][0] = 0;
    }
    else {
      v[i][0] = mid_prev;
      trace[i][0] = 1;
    }
    
    /* handle case of matching across
     * 3 cases: [\_,_|_]* OR [_|_,_|_]* OR [_/,_|_]
     */
    // CASE 1: [\_,_|_]
    double prev_mid = scores[i].mid + scores[i].skip + v[i - 1][0];
    // CASE 2: [_|_,_|_]*
    double mid_mid = scores[i].mid + scores[i].skip + v[i - 1][1];
    // CASE 3: [_/,_|_]*
    double next_mid = scores[i].mid + v[i - 1][2];
    
    if (prev_mid > mid_mid) {
      if (prev_mid > next_mid) {
	v[i][1] = prev_mid;
	trace[i][1] = 0;
	trace[i][3] = 1;
      }
      else {
	v[i][1] = next_mid;
	trace[i][1] = 2;
      }
    }
    else {
      if (mid_mid > next_mid) {
	v[i][1] = mid_mid;
	trace[i][1] = 1;
	trace[i][3] = 1;
      }
      else {
	v[i][1] = next_mid;
	trace[i][1] = 2;
      }
    }
    
    /* handle case of matching next
     * 3 cases: [\_,_/]* OR [_|_,_/]* OR [_/,_/]
     */
    // CASE 1: [\_,_/]*
    double prev_next = scores[i].next + scores[i].skip + v[i - 1][0];
    // CASE 2: [_|_,_/]*
    double mid_next = scores[i].next + scores[i].skip + v[i - 1][1];
    // CASE 3: [_/,_/]
    double next_next = scores[i].next + v[i - 1][2];
    
    if (prev_next > mid_next) {
      if (prev_next > next_next) {
	v[i][2] = prev_next;
	trace[i][2] = 0;
	trace[i][3] = 1;
      }
      else {
	v[i][2] = next_next;
	trace[i][2] = 2;
      }
    }
    else {
      if (mid_next > next_next) {
	v[i][2] = mid_next;
	trace[i][2] = 1;
	trace[i][3] = 1;
      }
      else {
	v[i][2] = next_next;
	trace[i][2] = 2;
      }
    }
  }

  // add the value of the final skip
  v.back()[0] += final_skip;
  v.back()[1] += final_skip;

  // If the final skip can be used, then use it
  if (max(v.back()[0], v.back()[1]) >= v.back()[2] && final_skip > 0) {
    path.push_back(make_pair(final_a, final_b));
    assert(path.back().first != Dat::NONE ||
	   path.back().second != Dat::NONE);
  }
  
  // Initialize the traceback
  size_t prev = 0;
  if (v.back()[0] > max(v.back()[1], v.back()[2])) {
    path.push_back(make_pair(scores.back().id, scores.back().prev_id));
    assert(path.back().first != Dat::NONE ||
	   path.back().second != Dat::NONE);
    prev = trace.back()[0];
  }
  else if (v.back()[1] > v.back()[2]) {
    path.push_back(make_pair(scores.back().id, scores.back().mid_id));
    assert(path.back().first != Dat::NONE ||
	   path.back().second != Dat::NONE);
    prev = trace.back()[1];
  }
  else {
    path.push_back(make_pair(scores.back().id, scores.back().next_id));
    assert(path.back().first != Dat::NONE ||
	   path.back().second != Dat::NONE);
    prev = trace.back()[2];
  }
  if (trace.back()[3]) {
    path.push_back(make_pair(scores.back().skip_id_a, scores.back().skip_id_b));
    assert(path.back().first != Dat::NONE ||
	   path.back().second != Dat::NONE);
  }  

  // Trace back, associating the pairs of segments
  for (size_t j = trace.size() - 1; j > 0; --j) {
    const size_t k = j - 1;
    if (prev == 0) {
      path.push_back(make_pair(scores[k].id, scores[k].prev_id));
      assert(path.back().first != Dat::NONE ||
	     path.back().second != Dat::NONE);
      prev = trace[k][0];
    }
    else if (prev == 1) {
      path.push_back(make_pair(scores[k].id, scores[k].mid_id));
      assert(path.back().first != Dat::NONE ||
	     path.back().second != Dat::NONE);
      prev = trace[k][1];
    }
    else {
      path.push_back(make_pair(scores[k].id, scores[k].next_id));
      assert(path.back().first != Dat::NONE ||
	     path.back().second != Dat::NONE);
      prev = trace[k][2];
    }
    if (trace[k][3]) {
      path.push_back(make_pair(scores[k].skip_id_a, 
			       scores[k].skip_id_b));
      assert(path.back().first != Dat::NONE ||
	     path.back().second != Dat::NONE);
    }
  }
  // Reverse the path, making sure it is in increasing order
  reverse(path.begin(), path.end());
}


void
elim_type(vector<Dat> &dats) {
  vector<Dat> t;
  size_t prev = 0;
  for (size_t i = 0; i < dats.size(); ++i) {
    if (!dats[i].type) {
      t.push_back(dats[i]);
      size_t best_idx = 0;
      double best_score = 0;
      for (size_t j = prev + 1; j < i; ++j)
	if (dats[j].mid > best_score) {
	  best_score = dats[j].mid;
	  best_idx = j;
	}
      if (best_score > 0) {
	t.back().skip = dats[best_idx].mid;
	t.back().skip_id_a = dats[best_idx].id;
	t.back().skip_id_b = dats[best_idx].mid_id;
      }
      prev = i;
    }
  }
  dats.swap(t);
}


void
partition_segments(const vector<SimpleGenomicRegion> &a,
		   const vector<SimpleGenomicRegion> &b,
		   vector<vector<SimpleGenomicRegion> > &a_out,
		   vector<vector<SimpleGenomicRegion> > &b_out) {
  
  vector<SimpleGenomicRegion> tmp(a);
  tmp.insert(tmp.end(), b.begin(), b.end());
  sort(tmp.begin(), tmp.end());
  
  collapse(tmp);

  a_out.resize(tmp.size());
  b_out.resize(tmp.size());
  size_t a_idx = 0, b_idx = 0;
  for (size_t i = 0; i < tmp.size(); ++i) {
    while (a_idx < a.size() && tmp[i].overlaps(a[a_idx])) {
      a_out[i].push_back(a[a_idx]);
      ++a_idx;
    }
    while (b_idx < b.size() && tmp[i].overlaps(b[b_idx])) {
      b_out[i].push_back(b[b_idx]);
      ++b_idx;
    }
  }
}


void
match_segments(const vector<SimpleGenomicRegion> &a_in,
	       const vector<SimpleGenomicRegion> &b_in,
// 	       vector<vector<GenomicRegion> > &a_parts,
// 	       vector<vector<GenomicRegion> > &b_parts,
	       vector<pair<size_t, size_t> > &matching) {
  
  vector<vector<SimpleGenomicRegion> > a;
  vector<vector<SimpleGenomicRegion> > b;
  partition_segments(a_in, b_in, a, b);
  
  size_t a_offset = 0;
  size_t b_offset = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    if (!a[i].empty() && !b[i].empty()) {
      
      // assign index to each genomic region
      vector<pair<size_t, SimpleGenomicRegion> > a_work;
      for (size_t j = 0; j < a[i].size(); ++j)
	a_work.push_back(make_pair(a_offset + j, a[i][j]));

      vector<pair<size_t, SimpleGenomicRegion> > b_work;
      for (size_t j = 0; j < b[i].size(); ++j)
	b_work.push_back(make_pair(b_offset + j, b[i][j]));
      
      // convert to scores
      vector<Dat> dats;
      get_triple_scores(a_work, b_work, dats);
      
      bool all_type = true;
      for (size_t j = 0; j < dats.size(); ++j)
	if (!dats[j].type)
	  all_type = false;
      
      if (!all_type) {
	
	Dat final(dats.back());
	if (!final.type)
	  final.mid = 0;
	
	elim_type(dats);
	
	// need to pass an offset for a and b to this function
	assign_parts(dats, final.mid, final.id, final.mid_id, matching);
      }
      else {
	
	size_t best_idx = 0;
	size_t best_other = 0;
	double best_score = 0;
	for (size_t j = 0; j < dats.size(); ++j)
	  if (dats[j].mid > best_score) {
	    best_score = dats[j].mid;
	    best_idx = dats[j].id;
	    best_other = dats[j].mid_id;
	  }
	
	matching.push_back(make_pair(best_idx, best_other));
	assert(matching.back().first != Dat::NONE ||
	       matching.back().second != Dat::NONE);
	
      }
    }
    a_offset += a[i].size();
    b_offset += b[i].size();
  }
  
  sort(matching.begin(), matching.end());
}

int
main(int argc, const char **argv) {
  
  try {
    
    /* FILES */
//     const char *bedfilea = 0, *bedfileb = 0;
//     const char *outfilea = 0, *outfileb = 0;
//     int VERBOSE = 0;
    
//     ////////////////////// COMMAND LINE OPTIONS /////////////////////////
//     static struct poptOption optionsTable[] = {
//       { "verbose", 'v', POPT_ARG_NONE, &VERBOSE, 0,
// 	"print more information" },
//       POPT_AUTOHELP POPT_TABLEEND
//     };
        
//     /****************** GET COMMAND LINE ARGUMENTS ***************************/
//     poptContext optCon = poptGetContext("bedmatch", argc, argv, optionsTable, 0);
//     poptSetOtherOptionHelp(optCon, "<regions_a> <regions_b> <outfile_a> <outfile_b>");
//     if (argc < 2) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_SUCCESS;
//     }
//     char c;
//     if ((c = poptGetNextOpt(optCon)) < -1) {
//       cerr << "bedmatch: bad argument "
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
//     if (poptPeekArg(optCon) == 0) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_FAILURE;
//     }
//     else outfilea = poptGetArg(optCon);
//     if (poptPeekArg(optCon) == 0) {
//       poptPrintHelp(optCon, stderr, 0);
//       return EXIT_FAILURE;
//     }
//     else outfileb = poptGetArg(optCon);
//     if (poptPeekArg(optCon) != 0) {
//       cerr << "bedmatch: leftover argument " << poptGetArg(optCon) << endl;
//       return EXIT_FAILURE;
//     }
//     poptFreeContext(optCon);

		  bool VERBOSE = false;
    OptionParser opt_parse("bedmatch", "A program for matching  two bed files",
			   "<File_A> <File_B> [Out_File_A] [Out_File_B]");
	opt_parse.add_opt("verbose", 'v', "print more information", false, VERBOSE );
	
	
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
	if (leftover_args.size() < 4)
	{
			cerr << "Need two file names" << endl;
			return EXIT_SUCCESS;
	}
	
    const string bedfilea = leftover_args[0];
	const string bedfileb = leftover_args[1];
	const string outfilea = leftover_args[2];
	const string outfileb = leftover_args[3];

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
    
    vector<pair<size_t, size_t> > matching;
    match_segments(regions_a, regions_b, matching);
    
    vector<SimpleGenomicRegion> matched_a, matched_b;

    double dist_to_ends = 0;
    double total_dist = 0;

    for (size_t i = 0; i < matching.size(); ++i) {
      if (matching[i].first != Dat::NONE && matching[i].second != Dat::NONE) {
	matched_a.push_back(regions_a[matching[i].first]);
	matched_b.push_back(regions_b[matching[i].second]);
	dist_to_ends += (max(matched_a.back().get_start(), matched_b.back().get_start()) - 
			 min(matched_a.back().get_start(), matched_b.back().get_start())) +
	  (max(matched_a.back().get_end(), matched_b.back().get_end()) - 
	   min(matched_a.back().get_end(), matched_b.back().get_end()));
	total_dist += matched_a.back().get_width() + matched_b.back().get_width();
      }
    }

    const double correction = max(regions_a.size(), regions_b.size())/
      static_cast<double>(matching.size());
    if (VERBOSE)
      cerr << "DIST_TO_ENDS=" << dist_to_ends << "\t" 
    	   << "TOTAL_DIST=" << total_dist << "\t" 
    	   << "ERROR=" << dist_to_ends/total_dist << "\t"
    	   << "CORRECTED_ERROR=" << correction*dist_to_ends/total_dist << endl;
    ofstream out(outfilea.c_str());
    copy(matched_a.begin(), matched_a.end(), 
    	 std::ostream_iterator<SimpleGenomicRegion>(out, "\n"));
    out.close();
    out.open(outfileb.c_str());
    copy(matched_b.begin(), matched_b.end(), 
    	 std::ostream_iterator<SimpleGenomicRegion>(out, "\n"));
    out.close();
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
