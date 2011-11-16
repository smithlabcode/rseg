/*
  Copyright (C) 2008 Cold Spring Harbor Laboratory
  Copyright (C) 2008 University of Southern California
  Authors: Andrew D. Smith, PhD
  
  This file is part of rmap.

  rmap is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  rmap is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with rmap; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef EPITYPE_HPP
#define EPITYPE_HPP

#include <vector>
#include <string>

#include "RNG.hpp"

class Epiread;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
///// EPITYPE CLASS

class Epitype {
public:
  Epitype(const std::string &tt) : t(tt) {}
  Epitype() {}
  
  std::string get_type() const {return t;}
  size_t get_length() const {return t.length();}

  bool operator<(const Epitype &rhs) const {return t < rhs.t;}
  bool operator==(const Epitype &rhs) const {return t == rhs.t;}
  std::string tostring() const;
  
  static void extract_epitypes(const std::vector<Epiread> &epireads,
			       std::vector<Epitype> &epitypes);
  bool is_methylated(size_t i) const {return t[i] == '1';}
  bool is_consistent(const std::string states, const size_t pos) const {
    return states == t.substr(pos, states.length());
  }
  
  double distance(const Epitype &rhs) const;
  
  static size_t
  count_consistent(const std::vector<Epitype> &epitypes,
		   const std::vector<Epiread> &epireads);
  
private:
  std::string t;
};

Epitype
random_epitype(const Runif &runif, const size_t l);

Epitype
random_epitype(const Runif &runif, const size_t l, const size_t m);

inline std::ostream&
operator<<(std::ostream &s, const Epitype &ept) {return s << ept.tostring();}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
///// EPIREAD CLASS

class Epiread {
public:
  Epiread() : p(std::numeric_limits<size_t>::max()), 
	      read_id(std::numeric_limits<size_t>::max()) {}
  Epiread(const std::string t_, const size_t p_, const size_t rid) : 
    t(t_), p(p_), read_id(rid) {}
  
  std::string tostring() const;
  bool fits(const Epitype &e) const;
  bool is_consistent(const Epitype &e) const;
  bool is_methylated(size_t i) const {return t[i] == '1';}
  
  size_t get_length() const {return t.length();}
  size_t get_read_id() const {return read_id;}
  size_t get_pos() const {return p;}
  std::string get_type() const {return t;}
  
  bool operator<(const Epiread &rhs) const {
    return ((p <  rhs.p) || 
	    (p == rhs.p && t.length() <  rhs.t.length()) ||
	    (p == rhs.p && t.length() == rhs.t.length() && t < rhs.t));}
  
private:
  std::string t;
  size_t p;
  size_t read_id;
};

inline std::ostream&
operator<<(std::ostream &s, const Epiread &epr) {return s << epr.tostring();}

void
sample_frags(const Runif &runif, const Epitype &e,
	     const size_t n_frags, const size_t min_frag_length,
	     const size_t frag_length,
	     std::vector<Epiread> &epireads);

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
///// EPITYPE DISTRO CLASS

class EpitypeDistro {
public:
  EpitypeDistro(size_t l) : 
    params(std::vector<std::pair<double, double> >(l, std::make_pair(0.5, 0.5))) {}
  EpitypeDistro(const Epitype &e, const double pseudo);
  void estimate_params(const std::vector<Epiread> &er, const std::vector<double> &probs);
  double log_likelihood(const Epiread &f) const;
  std::string tostring() const;
  Epitype to_epitype() const;
private:
  std::vector<std::pair<double, double> > params;

  static double min_prob;
};

inline std::ostream&
operator<<(std::ostream &s, const EpitypeDistro &epd) {return s << epd.tostring();}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
///// EPITYPE DISTRO EM CLASS

class EpitypeDistroEM {
public:
  
  static double
  expectation_maximization(const std::vector<Epiread> &epireads, 
			   std::vector<EpitypeDistro> &epi_distros,
			   std::vector<double> &mixing);
  
  static void
  set_parameters(const double t, const double mp, const size_t itr) {
    tolerance = t;
    min_prob = mp;
    max_iterations = itr;
  }
  
  static void 
  greedy_starting_point(const std::vector<Epitype> &all_epitypes, 
			const std::vector<Epiread> &epireads,
			const size_t n_true_types, 
			std::vector<EpitypeDistro> &epi_distros);
  
private:
  
  static double
  expectation_step(const std::vector<Epiread> &epireads, 
		   const std::vector<double> &mixing,
		   const std::vector<EpitypeDistro> &epi_distros, 
		   std::vector<std::vector<double> > &probs);
  
  static void
  maximization_step(const std::vector<Epiread> &epireads, 
		    const std::vector<std::vector<double> > &probs,
		    std::vector<double> &mixing, 
		    std::vector<EpitypeDistro> &epi_distros);
  
  static double tolerance;
  static double min_prob;
  static size_t max_iterations;
};

#endif
