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

#include "GenomicProfile.hpp"

#include <fstream>

using std::vector;
using std::string;


string
GenomicProfile::tostring() const {
  std::ostringstream ss;
  ss << "fixedStep"
     << " chrom=" << region.get_chrom() 
     << " start=" << region.get_start() + 1 // BECAUSE WIG FILES ARE 1 BASED!!!
     << " step=" << step << " span=" << step;
  for (vector<double>::const_iterator i(values.begin()); i != values.end(); ++i)
    ss << std::endl << *i;
  return ss.str();
}

std::ostream& 
operator<<(std::ostream& the_stream, const GenomicProfile& profile) {
  return the_stream << profile.tostring();
}

static void
parse_profile_name(const char *s, const size_t len,
		   string &chrom, size_t &start, size_t &step) {
  size_t i = 0;
  
  // the chrom
  while (s[i] != '=' && i < len) ++i;
  size_t j = i + 1;
  while (!isspace(s[i]) && i < len) ++i;
  chrom = string(s + j, i - j);
  
  // start of the region (a positive integer)
  while (s[i] != '=' && i < len) ++i;
  j = i + 1;
  while (!isspace(s[i]) && i < len) ++i;
  start = atoi(s + j);
  assert(start > 0);
  start -= 1; // BECAUSE WIG FILES ARE 1 BASED :-(

  // end of the region (a positive integer)
  while (s[i] != '=' && i < len) ++i;
  j = i + 1;
  while (!isspace(s[i]) && i < len) ++i;
  step = atoi(s + j);
}



static char *
find_next_end(char *c, char *buffer_end, char * &line_end, size_t &n_lines) {

  static const char *marker = "fixedStep";
  static const size_t lim = 9;

  // find the end of the current line
  line_end = std::find(c, buffer_end, '\n') + 1;
  c = line_end;
  assert(c < buffer_end);
  
  n_lines = 0;
  // loop until found another "\nfixedStep"
  while (c != buffer_end) {
    size_t j = 0;
    while (j < lim && c[j] == marker[j++]);
    // ++j;
    if (j == lim) return c;
    else {
      c = std::find(c, buffer_end, '\n') + 1;
      ++n_lines;
    }
  }
  return buffer_end;
}

static bool is_track_line(const char *c) {
  static const char *track = "track";
  static const size_t track_size = 5;
  for (size_t i = 0; i < track_size; ++i)
    if (c[i] != track[i]) return false;
  return true;
}

void
ReadWIGFile(string filename, vector<GenomicProfile> &prof) {
  
  std::ifstream in(filename.c_str());
  if (!in.good())
    throw WIGFileException("cannot open input file " + filename);
  
  size_t begin_pos = in.tellg();
  in.seekg(0, std::ios_base::end);
  size_t end_pos = in.tellg();
  in.seekg(0, std::ios_base::beg);
  
  size_t filesize = end_pos - begin_pos;
  char *buffer = new char[filesize + 1];
  
  in.read(buffer, filesize);
  prof.reserve(std::count(buffer, buffer + filesize, '\n'));

  string chrom;
  size_t start= 0;
  size_t step = 0;
  
  char *buffer_end = buffer + filesize;
  char *line_end = 0;
  size_t n_lines = 0;
  char *c = buffer;
  vector<double> vals;

  if (is_track_line(c))
    c = find_next_end(buffer, buffer_end, line_end, n_lines);
  
  while (c != buffer_end) {
    
    // Find the next end of a region
    char *end = find_next_end(c, buffer_end, line_end, n_lines);
    
    // parse the name line
    parse_profile_name(c, line_end - c, chrom, start, step);
    
    // parse the numbers
    c = line_end;
    vals.reserve(n_lines);
    while (c != end) {
      vals.push_back(atof(c));
      c = std::find(c, end, '\n') + 1;
    }
    prof.push_back(GenomicProfile(chrom, start, step, vals));
    vals.clear();
    c = end;
  }
  delete[] buffer;
  in.close();
}

void
WriteWIGFile(const string filename, const string track_name, 
	     const vector<GenomicProfile> &prof) {
  std::ofstream wigout(filename.c_str());
  wigout << "track type=wiggle_0";
  if (track_name.length() > 0)
    wigout << " name=" << track_name;
  wigout << std::endl;
  copy(prof.begin(), prof.end(), 
       std::ostream_iterator<GenomicProfile>(wigout, "\n"));
  wigout.close();
}

void
WriteWIGFile(const string filename, const string track_name, 
	     const GenomicProfile &prof) {
  std::ofstream wigout(filename.c_str());
  wigout << "track type=wiggle_0";
  if (track_name.length() > 0)
    wigout << " name=" << track_name << std::endl;
  wigout << std::endl;
  wigout << prof << std::endl;
  wigout.close();
}
