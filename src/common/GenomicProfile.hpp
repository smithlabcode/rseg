/*
 * Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song and Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#ifndef GENOMIC_PROFILE_HPP
#define GENOMIC_PROFILE_HPP

#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

class GenomicProfile {
public:
  
  typedef std::vector<double>::iterator pointer;
  typedef const std::vector<double>::const_iterator const_pointer;
  typedef double& reference;
  typedef const double & const_reference;
  typedef pointer iterator;
  typedef const_pointer const_iterator;
  
  GenomicProfile(const SimpleGenomicRegion &r,
		 const std::vector<double> &v, const size_t s) :
    region(r), values(v), step(s) {}
  GenomicProfile(const std::string &c, const size_t st, const size_t stp,
		 const std::vector<double> &v) :
    region(SimpleGenomicRegion(c, st, stp*v.size())), values(v), step(stp) {}
  GenomicProfile(const SimpleGenomicRegion &r, const size_t s) :
    region(r), step(s) {}
  std::string tostring() const;
  
  // accessors
  SimpleGenomicRegion get_region() const {return region;}
  std::vector<double> get_values() const {return values;}
  size_t get_step() const {return step;}
  size_t get_size() const {return values.size();}
  
  iterator begin() {return values.begin();}
  iterator end() {return values.end();}
  const_iterator begin() const {return values.begin();}
  const_iterator end() const {return values.end();}
  reference operator[](size_t n) {return *(values.begin() + n);}
  const_reference operator[](size_t n) const {return *(values.begin() + n);}
  
  // mutators
  void set_region(const SimpleGenomicRegion &r) {region = r;}
  void set_values(const std::vector<double> &v) {values = v;}
  void set_step(size_t s) {step = s;}
  
  void set_end(const size_t e) {region.set_end(e);}
  void set_end() {region.set_end(step*values.size());}
  
private:
  SimpleGenomicRegion region;
  std::vector<double> values;
  size_t step;
};

std::ostream& 
operator<<(std::ostream& the_stream, const GenomicProfile& profile);

std::ostream& 
operator<<(std::ostream& the_stream, const std::vector<GenomicProfile>& profile);

void
ReadWIGFile(std::string filename, std::vector<GenomicProfile> &the_profiles);
void
WriteWIGFile(const std::string filename, const std::string track_name, 
	     const std::vector<GenomicProfile> &prof);
void
WriteWIGFile(const std::string filename, const std::string track_name, 
	     const GenomicProfile &prof);

class GenomicProfileException : public SMITHLABException {
public:
  GenomicProfileException(std::string s = std::string()) : SMITHLABException(s) {}
};

class WIGFileException : public SMITHLABException {
public:
  WIGFileException(std::string s = std::string()) throw() : SMITHLABException(s) {}
};

#endif
