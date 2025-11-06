/* deadzones: A program for identifying genomic deadzones
 *
 * Copyright (C) 2009 University of Southern California and
 *                    Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <cmath>
#include <filesystem>
#include <numeric>
#include <unordered_map>

class IndexLess {
public:
  IndexLess(const std::size_t k, const std::string &s) :
    kmer(k), itr(s.begin()) {}
  bool
  operator()(std::size_t a, std::size_t b) const {
    const std::string::const_iterator lim(itr + a + kmer);
    std::string::const_iterator a_itr(itr + a), b_itr(itr + b);
    while (a_itr < lim && *(a_itr) == *(b_itr)) {
      ++a_itr;
      ++b_itr;
    }
    return (*a_itr < *b_itr && a_itr < lim);
  }

private:
  std::size_t kmer;
  const std::string::const_iterator itr;
};

template <class In>
bool
lexico_equal(In first, In last, In first2) {
  while (first != last)
    if (*first++ != *first2++)
      return false;
  return true;
}

static void
sort_index(const bool VERBOSE, const std::size_t kmer,
           const std::string &prefix, const std::string &seq,
           std::vector<std::size_t> &ambigs,
           const std::unordered_map<std::size_t, std::size_t> &invalid_pool) {

  if (VERBOSE)
    std::cerr << "[BUILDING INDEX] ";
  std::vector<std::size_t> index;
  const std::string::const_iterator lim(seq.end() - kmer + 1);
  for (std::string::const_iterator j = seq.begin(); j != lim; ++j)
    if ((lexico_equal(prefix.begin(), prefix.end(), j)) &&
        (!(invalid_pool.find(j - seq.begin()) != invalid_pool.end())))
      index.push_back(j - seq.begin());

  if (!index.empty()) {

    if (VERBOSE)
      std::cerr << "[SORTING INDEX] ";
    IndexLess index_less(kmer, seq);
    sort(index.begin(), index.end(), index_less);

    if (VERBOSE)
      std::cerr << "[FINDING DEADS] ";
    const std::size_t len = seq.length();
    const std::string::const_iterator start(seq.begin());
    const std::string::const_iterator end(start + kmer);
    std::size_t prev = index.front();
    bool prev_inserted = false;
    for (std::size_t i = 1; i < index.size(); ++i) {
      const std::size_t curr = index[i];
      if (lexico_equal(start + prev, end + prev, start + curr) &&
          prev + curr != len) {
        if (!prev_inserted)
          ambigs.push_back(prev);
        ambigs.push_back(curr);
        prev_inserted = true;
      }
      else
        prev_inserted = false;
      prev = curr;
    }
  }
  else if (VERBOSE)
    std::cerr << "[EMPTY INDEX] ";
}

static void
sort_index(const bool VERBOSE, const bool BISULFITE, const bool AG_WILDCARD,
           const std::size_t kmer, const std::size_t prefix_len,
           const std::string &seq, std::vector<std::size_t> &ambigs,
           const std::unordered_map<std::size_t, std::size_t> &invalid_pool) {

  static const float DENOM = CLOCKS_PER_SEC;

  const std::size_t n_prefix =
    static_cast<std::size_t>(pow(smithlab::alphabet_size, prefix_len));
  for (std::size_t i = 0; i < n_prefix; ++i) {
    const std::string prefix(i2mer(prefix_len, i));
    if (!BISULFITE ||
        ((!AG_WILDCARD && prefix.find('C') == std::string::npos) ||
         (AG_WILDCARD && prefix.find('G') == std::string::npos))) {
      const clock_t start(clock());
      if (VERBOSE)
        std::cerr << "[PREFIX=" << prefix << "] ";
      sort_index(VERBOSE, kmer, prefix, seq, ambigs, invalid_pool);
      const clock_t end(clock());
      if (VERBOSE)
        std::cerr << "[" << (end - start) / DENOM << " SEC] [DONE]" << '\n';
    }
  }
}

static void
write_dead(std::ofstream &out, const std::string &chrom_name, const char strand,
           std::vector<std::size_t>::const_iterator curr,
           const std::vector<std::size_t>::const_iterator lim) {
  assert(curr <= lim);
  std::size_t prev_ambig = *curr;
  ++curr;
  for (; curr < lim; ++curr)
    if (*curr - 1 != *(curr - 1)) {
      out << GenomicRegion(chrom_name, prev_ambig, *(curr - 1) + 1, "X", 0,
                           strand)
          << '\n';
      prev_ambig = *curr;
    }
  out << GenomicRegion(chrom_name, prev_ambig, *(curr - 1) + 1, "X", 0, strand)
      << '\n';
}

static void
write_dead(std::ofstream &out, const std::string &chrom_name, const char strand,
           std::vector<std::size_t>::const_iterator curr,
           const std::vector<std::size_t>::const_iterator lim,
           const std::vector<SimpleGenomicRegion> &gaps) {
  assert(curr <= lim);
  assert(gaps.size() == 0 || gaps.front().get_chrom() == chrom_name);

  std::vector<GenomicRegion> deadzones;
  std::size_t prev_ambig = *curr;
  ++curr;
  for (; curr < lim; ++curr)
    if (*curr - 1 != *(curr - 1)) {
      deadzones.push_back(
        GenomicRegion(chrom_name, prev_ambig, *(curr - 1) + 1, "X", 0, strand));
      prev_ambig = *curr;
    }
  deadzones.push_back(
    GenomicRegion(chrom_name, prev_ambig, *(curr - 1) + 1, "X", 0, strand));

  std::size_t i = 0, j = 0;
  while (i < deadzones.size() && j < gaps.size()) {
    const std::size_t s =
      std::max(deadzones[i].get_start(), gaps[j].get_start());
    const std::size_t e = std::min(deadzones[i].get_end(), gaps[j].get_end());

    if (s > e) {  //  no overlapping
      if (deadzones[i].get_end() < gaps[j].get_start()) {
        out << deadzones[i] << '\n';
        ++i;
      }
      else {
        out << gaps[j] << "\tT" << "\t" << 0 << "\t" << strand << '\n';
        ++j;
      }
    }
    else {
      deadzones[i].set_start(
        std::min(deadzones[i].get_start(), gaps[j].get_start()));
      deadzones[i].set_end(std::max(deadzones[i].get_end(), gaps[j].get_end()));
      ++j;
      while (i + 1 < deadzones.size() &&
             deadzones[i].get_end() >= deadzones[i + 1].get_start()) {
        ++i;
        deadzones[i].set_start(deadzones[i - 1].get_start());
        deadzones[i].set_end(
          std::max(deadzones[i].get_end(), deadzones[i - 1].get_end()));
      }
    }
  }
  while (i < deadzones.size()) {
    out << deadzones[i] << '\n';
    ++i;
  }
  while (j < gaps.size()) {
    out << gaps[j] << "\tT" << "\t" << 0 << "\t" << strand << '\n';
    ++j;
  }
}

static void
get_dead(const std::string &outfile, const std::size_t kmer,
         const std::vector<std::size_t> &seqoffsets,
         const std::vector<std::string> &chrom_names,
         std::vector<std::size_t> &ambigs,
         std::vector<std::vector<SimpleGenomicRegion>> &gaps) {

  const std::size_t max_offset = seqoffsets.back();
  for (std::size_t i = 0; i < ambigs.size(); ++i) {
    if (ambigs[i] >= max_offset)
      ambigs[i] = 2 * max_offset - ambigs[i] - kmer;
    assert(ambigs[i] < max_offset);
  }
  sort(ambigs.begin(), ambigs.end());
  ambigs.erase(std::unique(ambigs.begin(), ambigs.end()), ambigs.end());

  std::vector<std::size_t> offset_idx;
  std::size_t n_ambigs = ambigs.size();
  for (std::size_t i = 0, j = 0; i < seqoffsets.size() && j < n_ambigs; ++i) {
    while (j < n_ambigs && ambigs[j] < seqoffsets[i])
      ++j;
    offset_idx.push_back(j);
  }

  std::size_t total_length = 0;
  n_ambigs = ambigs.size();
  for (std::size_t i = 0, prev_idx = 0; i < offset_idx.size(); ++i) {
    for (std::size_t j = prev_idx; j < offset_idx[i]; ++j) {
      assert(j < n_ambigs);
      ambigs[j] -= total_length;
    }
    prev_idx = offset_idx[i];
    total_length = seqoffsets[i];
  }

  std::ofstream out(outfile);
  for (std::size_t i = 0, prev_idx = 0; i < offset_idx.size(); ++i) {
    write_dead(out, chrom_names[i], '+', ambigs.begin() + prev_idx,
               ambigs.begin() + offset_idx[i], gaps[i]);
    prev_idx = offset_idx[i];
  }
  out.close();
}

static void
get_dead_bs(const bool VERBOSE, const std::string &outfile,
            const std::size_t kmer, const std::vector<std::size_t> &seqoffsets,
            const std::vector<std::string> &chrom_names,
            std::vector<std::size_t> &ambigs) {
  assert(!ambigs.empty());
  sort(ambigs.begin(), ambigs.end());

  const std::size_t max_offset = seqoffsets.back();
  if (VERBOSE)
    std::cerr << "[PREPARING POS-STRAND BS DEADS]" << '\n';

  // Do the positive strand bisulfite deadzones
  const std::size_t lim =
    lower_bound(ambigs.begin(), ambigs.end(), max_offset) - ambigs.begin();

  // make a partition vector of the offsets, the last being "lim"
  std::vector<std::size_t> offset_idx;
  std::size_t n_ambigs = ambigs.size();
  for (std::size_t i = 0, j = 0; i < seqoffsets.size() && j < n_ambigs; ++i) {
    while (j < n_ambigs && ambigs[j] < seqoffsets[i])
      ++j;
    offset_idx.push_back(j);
  }

  std::size_t total_length = 0;
  for (std::size_t i = 0, prev_idx = 0; i < offset_idx.size(); ++i) {
    for (std::size_t j = prev_idx; j < offset_idx[i]; ++j)
      ambigs[j] -= total_length;
    prev_idx = offset_idx[i];
    total_length = seqoffsets[i];
  }

  std::ofstream out(outfile);

  for (std::size_t i = 0, prev_idx = 0; i < offset_idx.size(); ++i) {
    write_dead(out, chrom_names[i], '+', ambigs.begin() + prev_idx,
               ambigs.begin() + offset_idx[i]);
    prev_idx = offset_idx[i];
  }

  if (VERBOSE)
    std::cerr << "[PREPARING NEG-STRAND BS DEADS]" << '\n';
  // Move the negative strand deadzones into the first portion of the
  // vector and correct their indexes.
  for (std::size_t j = 0, i = lim; i < ambigs.size(); ++i)
    ambigs[j++] = 2 * max_offset - ambigs[i] - kmer;
  ambigs.erase(ambigs.end() - lim, ambigs.end());
  reverse(ambigs.begin(), ambigs.end());

  offset_idx.clear();
  n_ambigs = ambigs.size();
  for (std::size_t i = 0, j = 0; i < seqoffsets.size() && j < n_ambigs; ++i) {
    while (j < n_ambigs && ambigs[j] < seqoffsets[i])
      ++j;
    offset_idx.push_back(j);
  }

  total_length = 0;
  for (std::size_t i = 0, prev_idx = 0; i < offset_idx.size(); ++i) {
    for (std::size_t j = prev_idx; j < offset_idx[i]; ++j)
      ambigs[j] -= total_length;
    prev_idx = offset_idx[i];
    total_length = seqoffsets[i];
  }

  for (std::size_t i = 0, prev_idx = 0; i < offset_idx.size(); ++i) {
    write_dead(out, chrom_names[i], '-', ambigs.begin() + prev_idx,
               ambigs.begin() + offset_idx[i]);
    prev_idx = offset_idx[i];
  }
  out.close();
}

// This function appends the reverse complement in a space efficient way
static void
append_revcomp(std::string &long_seq) {
  const std::size_t seqlen = long_seq.length();
  long_seq.resize(2 * seqlen);
  copy(long_seq.begin(), long_seq.begin() + seqlen, long_seq.begin() + seqlen);
  revcomp_inplace(long_seq.begin() + seqlen, long_seq.end());
}

static void
identify_chromosomes(const bool VERBOSE, const std::string fasta_suffix,
                     const std::string chrom_file,
                     std::vector<std::string> &chrom_files) {
  if (VERBOSE)
    std::cerr << "[IDENTIFYING CHROMS] ";
  if (std::filesystem::is_directory(chrom_file))
    read_dir(chrom_file, fasta_suffix, chrom_files);
  else
    chrom_files.push_back(chrom_file);
  std::sort(chrom_files.begin(), chrom_files.end());
  if (VERBOSE) {
    std::cerr << "[DONE]" << '\n'
              << "chromosome files found (approx size):" << '\n';
    for (std::vector<std::string>::const_iterator i = chrom_files.begin();
         i != chrom_files.end(); ++i)
      std::cerr << *i << " (" << roundf(get_filesize(*i) / 1e06) << "Mbp)"
                << '\n';
    std::cerr << '\n';
  }
}

int
main(int argc, char *argv[]) {

  try {

    // Parameter variables
    std::size_t kmer = 0;
    std::size_t prefix_len = 5;
    std::string outfile;
    std::string fasta_suffix = "fa";

    bool VERBOSE = false;
    bool BISULFITE = false;
    bool AG_WILDCARD = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("deadzones", "program for finding deadzones",
                           "<1-or-more-FASTA-chrom-files>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      true, outfile);
    opt_parse.add_opt("kmer", 'k', "Width of k-mers", true, kmer);
    opt_parse.add_opt("prefix", 'p', "prefix length (default 5)", false,
                      prefix_len);
    // opt_parse.add_opt("bisulfite", 'B', "get bisulfite deadzones",
    //        false, BISULFITE);
    // opt_parse.add_opt("ag-wild", 'A', "A/G wildcard for bisulfite",
    //        false, AG_WILDCARD);
    opt_parse.add_opt("suffix", 's',
                      "suffix of FASTA files "
                      "(assumes -c indicates dir)",
                      false, fasta_suffix);
    opt_parse.add_opt("verbose", 'v', "print more run information", false,
                      VERBOSE);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string chrom_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::vector<std::string> seqfiles;
    identify_chromosomes(VERBOSE, fasta_suffix, chrom_file, seqfiles);

    std::string long_seq;
    std::vector<std::size_t> seqoffsets;
    std::vector<std::string> chrom_names;

    if (VERBOSE)
      std::cerr << "[READING SEQUENCE FILES]" << '\n';
    std::vector<std::vector<SimpleGenomicRegion>> gaps;
    for (std::size_t i = 0; i < seqfiles.size(); ++i) {
      if (std::filesystem::is_directory(seqfiles[i]))
        throw std::runtime_error("\"" + seqfiles[i] +
                                 "\" not a FASTA format sequence file?");
      std::vector<std::string> names, sequences;
      read_fasta_file(seqfiles[i], names, sequences);
      for (std::size_t j = 0; j < sequences.size(); ++j) {
        gaps.push_back(std::vector<SimpleGenomicRegion>());
        transform(sequences[j].begin(), sequences[j].end(),
                  sequences[j].begin(),
                  [&](const char x) { return std::toupper(x); });
        std::size_t start = sequences[j].find_first_of('N');
        while (start != std::string::npos) {
          const static std::size_t min_gap_size = 50;
          std::size_t end = sequences[j].find_first_not_of('N', start);
          end = end == std::string::npos ? sequences[j].size() : end;
          if (start + min_gap_size <= end)
            gaps.back().push_back(SimpleGenomicRegion(names[j], start, end));
          start = sequences[j].find_first_of('N', end);
        }
        long_seq += sequences[j];
        seqoffsets.push_back(long_seq.length());
        chrom_names.push_back(names[j]);
      }
      if (VERBOSE)
        std::cerr << seqfiles[i] << "\t(SEQS: " << std::size(names) << ")\n";
    }

    if (VERBOSE)
      std::cerr << "[PREPARING CONCATENATED SEQUENCE]\n";
    append_revcomp(long_seq);

    if (BISULFITE) {
      if (AG_WILDCARD)
        std::replace(std::begin(long_seq), std::end(long_seq), 'G', 'A');
      else
        std::replace(std::begin(long_seq), std::end(long_seq), 'C', 'T');
    }

    if (VERBOSE)
      std::cerr << "[PREPARING INVALID INDEXES]\n";
    std::unordered_map<std::size_t, std::size_t> invalid_pool;
    std::size_t max = seqoffsets[seqoffsets.size() - 1];
    for (std::size_t i = 0; i < seqoffsets.size(); i++) {
      for (std::size_t j = seqoffsets[i] - kmer + 1; j <= seqoffsets[i] - 1;
           j++)
        invalid_pool[j] = 1;
      for (std::size_t j = max + (max - seqoffsets[i]) - kmer + 1;
           j <= max + (max - seqoffsets[i] - 1); j++)
        invalid_pool[j] = 1;
    }

    if (VERBOSE)
      std::cerr << "[IDENTIFYING AMBIGUOUS INDEXES]\n";
    std::vector<std::size_t> ambigs;

    sort_index(VERBOSE, BISULFITE, AG_WILDCARD, kmer, prefix_len, long_seq,
               ambigs, invalid_pool);

    long_seq.clear();

    if (ambigs.empty()) {
      if (VERBOSE)
        std::cerr << "[NO DEADZONES FOUND]\n";
    }
    else {
      if (BISULFITE) {
        if (VERBOSE)
          std::cerr << "[PREPARING BS DEADZONES]\n";
        get_dead_bs(VERBOSE, outfile, kmer, seqoffsets, chrom_names, ambigs);
      }
      else {
        if (VERBOSE)
          std::cerr << "[PREPARING DEADZONES]" << '\n';
        get_dead(outfile, kmer, seqoffsets, chrom_names, ambigs, gaps);
      }
    }
  }
  catch (std::runtime_error &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
