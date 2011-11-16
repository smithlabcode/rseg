#include "Epitype.hpp"

#include <sstream>
#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::cerr;
using std::endl;
using std::sort;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

class EpiTreeNode {
public:
  EpiTreeNode() : weight(0), zero(0), one(0) {}
  ~EpiTreeNode();
  string tostring(size_t depth) const;
  void extract_epitypes(const string ep_prefix, 
			const double prefix_score,
			vector<Epitype> &epitypes) const;
private:
  void insert(const string &s, size_t pos, size_t depth, bool first = false);
  double weight;
  EpiTreeNode *zero;
  EpiTreeNode *one;
  
  static bool is_zero(char c) {return c == '0';}
  static bool is_one(char c) {return c == '1';}
  static bool node_exists(const EpiTreeNode *n) {return n;}
  
  friend class EpiTree;
};

EpiTreeNode::~EpiTreeNode() {
  if (zero) {
    delete zero;
    zero = 0;
  }
  if (one) {
    delete one;
    one = 0;
  }
}

void
EpiTreeNode::extract_epitypes(const string ep_prefix, const double prefix_score,
			      vector<Epitype> &epitypes) const {
  if (!node_exists(zero) && !node_exists(one))
    epitypes.push_back(Epitype(ep_prefix));
  if (node_exists(zero)) 
    zero->extract_epitypes(ep_prefix + "0", prefix_score + zero->weight, epitypes);
  if (node_exists(one)) 
    one->extract_epitypes(ep_prefix + "1", prefix_score + one->weight, epitypes);
}

string
EpiTreeNode::tostring(size_t depth) const {
  std::ostringstream ss;
  if (node_exists(zero))
    ss << string(depth, ' ') << "0[" 
       << zero->weight << "]" << std::endl
       << zero->tostring(depth + 1);
  if (node_exists(one))
    ss << string(depth, ' ') 
       << "1[" << one->weight << "]" << std::endl
       << one->tostring(depth + 1);
  return ss.str();
}

std::ostream&
operator<<(std::ostream &ss, const EpiTreeNode &n) {return ss << n.tostring(0);}

void
EpiTreeNode::insert(const string &s, size_t pos, size_t depth, bool first) {
  if (depth > 0) {
    if (node_exists(zero))
      zero->insert(s, 0, depth - 1);
    if (node_exists(one))
      one->insert(s, 0, depth - 1);
  }
  else {
    weight += 1;
    if (is_zero(s[pos])) {
      if (!node_exists(zero) && (first || pos > 0))
	zero = new EpiTreeNode();
      if (node_exists(zero))
	zero->insert(s, pos + 1, 0);
    }
    if (is_one(s[pos])) {
      if (!node_exists(one) && (first || pos > 0))
	one = new EpiTreeNode();
      if (node_exists(one))
	one->insert(s, pos + 1, 0);
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

class EpiTree {
public:
  EpiTree() : root(new EpiTreeNode) {}
  string tostring() const;
  void extract_epitypes(vector<Epitype> &epitypes) const;
  void insert(const string &s, size_t depth);
private:
  EpiTreeNode *root;
};

void
EpiTree::insert(const string &s, size_t depth) {
  if (!EpiTreeNode::node_exists(root)) root = new EpiTreeNode;
  root->insert(s, 0, depth, true);
}

void
EpiTree::extract_epitypes(vector<Epitype> &epitypes) const {
  if (EpiTreeNode::node_exists(root))
    root->extract_epitypes("", 0, epitypes);
  size_t max_len = 0;
  for (size_t i = 0; i < epitypes.size(); ++i)
    max_len = std::max(max_len, epitypes[i].get_length());
  size_t j = 0;
  for (size_t i = 0; i < epitypes.size(); ++i)
    if (epitypes[i].get_length() == max_len)
      epitypes[j++] = epitypes[i];
  epitypes.erase(epitypes.begin() + j, epitypes.end());
}

string
EpiTree::tostring() const {
  return (EpiTreeNode::node_exists(root)) ? root->tostring(0) : string();
}
std::ostream&
operator<<(std::ostream &ss, const EpiTree &t) {return ss << t.tostring();}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



string
Epitype::tostring() const {
//   std::ostringstream ss;
//   ss << t;
  return t;//ss.str();
}



double
Epitype::distance(const Epitype &rhs) const {
  size_t diffs = 0;
  for (size_t i = 0; i < t.length(); ++i)
    diffs += (rhs.t[i] != t[i]);
  return diffs;
}




void
Epitype::extract_epitypes(const vector<Epiread> &epireads,
			  vector<Epitype> &epitypes) {
  using std::pair;
  vector<pair<size_t, string> > sorter;
  for (size_t i = 0; i < epireads.size(); ++i)
    sorter.push_back(make_pair(epireads[i].get_pos(), epireads[i].get_type()));
  sort(sorter.begin(), sorter.end());

  EpiTree tree;
  for (size_t i = 0; i < sorter.size(); ++i)
    tree.insert(sorter[i].second, sorter[i].first);
  tree.extract_epitypes(epitypes);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

Epitype
random_epitype(const Runif &runif, const size_t l) {
  string s;
  for (size_t i = 0; i < l; ++i)
    s += (runif.runif(0.0,1.0) > 0.5) ? '1' : '0';
  return Epitype(s);
}


Epitype
random_epitype(const Runif &runif, const size_t l, const size_t m) {
  string s;
  const size_t limit = runif.runif(l, m);
  for (size_t i = 0; i < limit; ++i)
    s += (runif.runif(0.0,1.0) > 0.5) ? '1' : '0';
  return Epitype(s);
}


void
sample_frags(const Runif &runif, const Epitype &e,
	     const size_t n_frags, const size_t min_frag_length,
	     const size_t max_frag_length,
	     vector<Epiread> &epireads) {
  vector<pair<size_t, string> > sorter;
  const string s(e.get_type());
  const size_t lim = s.length();
  for (size_t i = 0; i < n_frags; ++i) {
    const size_t current_frag_length = runif.runif(min_frag_length, max_frag_length + 1);
    cerr << current_frag_length << endl;
    const size_t the_frag = runif.runif(0ul, lim - current_frag_length + 1);
    sorter.push_back(make_pair(the_frag, s.substr(the_frag, current_frag_length)));
  }
  sort(sorter.begin(), sorter.end());
  for (size_t i = 0; i < sorter.size(); ++i)
    epireads.push_back(Epiread(sorter[i].second, sorter[i].first, 0));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

string
Epiread::tostring() const {
  std::ostringstream ss;
  ss << t << "\t" << p << "\t" << read_id;
  return ss.str();
}

bool
Epiread::fits(const Epitype &e) const {
  return p + t.length() <= e.get_length();
}

bool
Epiread::is_consistent(const Epitype &e) const {
  return fits(e) && e.is_consistent(t, p);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

double EpitypeDistro::min_prob = 1e-6;

void
EpitypeDistro::estimate_params(const vector<Epiread> &er, const vector<double> &probs) {
  for (size_t i = 0; i < er.size(); ++i) {
    assert(i < probs.size() && probs[i] >= 0 && probs[i] <= 1);
    
    const size_t offset = er[i].get_pos();
    const size_t epiread_len = er[i].get_length();
    for (size_t j = 0; j < epiread_len; ++j) {
      if (er[i].is_methylated(j))
	params[offset + j].first += probs[i];
      else params[offset + j].second += probs[i];
    }
  }
  for (size_t i = 0; i < params.size(); ++i) {
    const double denom = params[i].first + params[i].second;
    params[i].first /= denom;
    params[i].second /= denom;
  }
}

double
EpitypeDistro::log_likelihood(const Epiread &e) const {
  const size_t k = e.get_pos();
  const size_t epiread_len = e.get_length();
  double score = 0;
  for (size_t i = 0; i < epiread_len; ++i) {
    score += std::log((e.is_methylated(i) ? params[k + i].first : params[k + i].second));
    assert(finite(score));
  }
  return score;
}

#include <iomanip>

string
EpitypeDistro::tostring() const {
  std::ostringstream ss;
  if (!params.empty())
    ss << std::setprecision(3) << std::fixed
       << params.front().first << "\t" << params.front().second;
  for (size_t i = 1; i < params.size(); ++i)
    ss << std::endl << params[i].first << "\t" << params[i].second;
  return ss.str();
}

EpitypeDistro::EpitypeDistro(const Epitype &e, const double pseudo) {
  const size_t epitype_len = e.get_length();
  params.resize(epitype_len, make_pair(pseudo, pseudo));
  for (size_t i = 0; i < epitype_len; ++i) {
    if (e.is_methylated(i)) params[i].first += 1;
    else params[i].second += 1;
    const double denom = params[i].first + params[i].second;
    params[i].first /= denom;
    params[i].second /= denom;
  }
  for (size_t i = 0; i < params.size(); ++i) {
    assert(params[i].first >= 0 && params[i].first <= 1);
    assert(params[i].second >= 0 && params[i].second <= 1);
  }
}

Epitype
EpitypeDistro::to_epitype() const {
  string s(params.size(), '0');
  for (size_t i = 0; i < params.size(); ++i)
    if (params[i].first > 0.5) s[i] = '1';
  return Epitype(s);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

double EpitypeDistroEM::tolerance = 1e-10;
double EpitypeDistroEM::min_prob = 1e-10;
size_t EpitypeDistroEM::max_iterations = 100;


static double
epitype_log_sum_log_vec(const vector<double> &vals, size_t limit) {
  const vector<double>::const_iterator x = 
    max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i)
    if (i != max_idx)
      sum += exp(vals[i] - max_val);
  return max_val + std::log(sum);
}

static double
epitype_log_sum_log_vec(const vector<vector<double> > &vals, size_t limit, size_t index) {
  double max_val = -std::numeric_limits<double>::max();
  size_t max_idx = 0;
  for (size_t i = 0; i < limit; ++i) {
    if (vals[i][index] > max_val) {
      max_val = vals[i][index];
      max_idx = i;
    }
  }
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i)
    if (i != max_idx)
      sum += exp(vals[i][index] - max_val);
  return max_val + std::log(sum);
}

double
EpitypeDistroEM::expectation_step(const vector<Epiread> &epireads, 
				  const vector<double> &mixing,
				  const vector<EpitypeDistro> &epi_distros, 
				  vector<vector<double> > &probs) {
  vector<double> log_mixing(mixing);
  transform(log_mixing.begin(), log_mixing.end(), log_mixing.begin(),
	    std::ptr_fun<double, double>(&std::log));
  
  const size_t n_reads = epireads.size();
  const size_t n_epitypes = epi_distros.size();
  
  double score = 0;
  for (size_t i = 0; i < n_reads; ++i) {
    for (size_t j = 0; j < n_epitypes; ++j) {
      const double current_part = log_mixing[j] + epi_distros[j].log_likelihood(epireads[i]);
      assert(finite(current_part));
      probs[j][i] = std::max(current_part, log(min_prob));
    }
    const double denom = epitype_log_sum_log_vec(probs, probs.size(), i);
    for (size_t j = 0; j < n_epitypes; ++j)
      probs[j][i] = exp(probs[j][i] - denom);
    score += denom;
  }
  return score;
}

void
EpitypeDistroEM::maximization_step(const vector<Epiread> &epireads, 
				   const vector<vector<double> > &probs,
				   vector<double> &mixing, 
				   vector<EpitypeDistro> &epi_distros) {
  const size_t n_reads = epireads.size();
  const size_t n_epitypes = epi_distros.size();
  
  for (size_t i = 0; i < n_epitypes; ++i)
    epi_distros[i].estimate_params(epireads, probs[i]);
  
  vector<double> log_mixing(n_epitypes);
  for (size_t i = 0; i < n_epitypes; ++i)
    log_mixing[i] = epitype_log_sum_log_vec(probs[i], n_reads);
  
  const double mix_sum = epitype_log_sum_log_vec(log_mixing, n_epitypes);
  for (size_t i = 0; i < n_epitypes; ++i)
    mixing[i] = exp(log_mixing[i] - mix_sum);
}

double
EpitypeDistroEM::expectation_maximization(const vector<Epiread> &epireads, 
					  vector<EpitypeDistro> &epi_distros,
					  vector<double> &mixing) {
  
  const size_t n_reads = epireads.size();
  const size_t n_distros = epi_distros.size();
  
  mixing.resize(n_distros, 1.0/n_distros);
  vector<vector<double> > probs(n_distros, vector<double>(n_reads, 0));
  
  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    const double score = expectation_step(epireads, mixing, epi_distros, probs);
    maximization_step(epireads, probs, mixing, epi_distros);
    if ((prev_score - score)/prev_score < tolerance)
      break;
    prev_score = score;
  }
  return prev_score;
}


static void
find_best_cover(const vector<Epitype> &epitypes, const vector<Epiread> &epireads,
		const vector<bool> &current_used, Epitype &best_epitype,
		vector<bool> &best_used) {
  
  const size_t n_reads = epireads.size();
  const size_t n_epitypes = epitypes.size();
  
  double best_score = 0;
  size_t best_idx = 0;
  for (size_t i = 0; i < n_epitypes; ++i) {
    double score = 0;
    for (size_t j = 0; j < n_reads; ++j)
      score += (!current_used[j] && epireads[j].is_consistent(epitypes[i]));
    if (score > best_score) {
      best_score = score;
      best_idx = i;
    }
  }
  best_epitype = epitypes[best_idx];
  
  best_used.resize(n_reads);
  for (size_t i = 0; i < n_reads; ++i)
    best_used[i] = (current_used[i] || epireads[i].is_consistent(best_epitype)); 
}

void
EpitypeDistroEM::greedy_starting_point(const vector<Epitype> &all_epitypes, 
				       const vector<Epiread> &epireads,
				       const size_t n_true_types, 
				       vector<EpitypeDistro> &epi_distros) {
  vector<bool> used(epireads.size(), false);
  for (size_t i = 0; i < n_true_types; ++i) {
    Epitype best_epitype;
    find_best_cover(all_epitypes, epireads, used, best_epitype, used);
    epi_distros.push_back(EpitypeDistro(best_epitype, min_prob));
  }
}

size_t
Epitype::count_consistent(const std::vector<Epitype> &epitypes,
			  const std::vector<Epiread> &epireads) {
  
  size_t consistent = 0;
  for (size_t i = 0; i < epireads.size(); ++i) {
    size_t j = 0;
    while (j < epitypes.size() && !epireads[i].is_consistent(epitypes[j]))
      ++j;
    if (j < epitypes.size())
      ++consistent;
  }
  return consistent;
}
