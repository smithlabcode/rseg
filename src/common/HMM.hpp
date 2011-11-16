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

#ifndef HMM_HEADER
#define HMM_HEADER

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <utility>
#include <iterator>

#include <cmath>
#include <cassert>

#include <gsl/gsl_sf_log.h>

#include "OptionParser.hpp"
#include "Distro.hpp"
#include "SplitDistro.hpp"
#include "smithlab_utils.hpp"
#include "RNG.hpp"

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::cin;
using std::getline;
using std::ostream_iterator;
using std::copy;
using std::pair;


class HMM
{
public:
    // constructors
    HMM() {}
    HMM(const size_t n);

    // simulators
    void simulate(const size_t state, pair<double, double> &obs) const;
    void simulate(const size_t n, vector<size_t> &states,
                  vector< pair<double, double> > &obs) const;
    void simulate(const vector<size_t> &states,
                  vector< pair<double, double> > &obs) const;

    // compute likelihood given sequence of observations
    double get_log_likelihood(const size_t &state, const  pair<double, double>  &obs) const;
    double	get_log_likelihood(const size_t &state,
                               const vector< pair<double, double> > &obs) const;
    double get_log_likelihood(const vector<size_t> &states,
                              const vector< pair<double, double> > &obs) const;

    // accessors: expensive functions!
    size_t get_states_n() const {return n_states;}
    vector<SplitDistro> get_emission_distros() const {return emis_distros;}
    vector<double> get_start_trans() const {return start_trans;}
    vector< vector<double> > get_trans() const  {return trans;}
    vector<double> get_end_trans() const {return end_trans;}

    // mutators
    void set_states_n(const size_t n);
    void set_emission_distros(const vector<SplitDistro> &distros)
    {
        copy(distros.begin(), distros.end(), emis_distros.begin());
    }

    void set_start_trans(const vector<double> &start_trans_p)
    {
        copy(start_trans_p.begin(), start_trans_p.end(),
             start_trans.begin());
    }

    void set_trans(const vector< vector<double> > &trans_p)
    {
        copy(trans_p.begin(), trans_p.end(),
             trans.begin());
    }

    void set_end_trans(const vector<double> &end_trans_p)
    {
        copy(end_trans_p.begin(), end_trans_p.end(),
             end_trans.begin());
    }

private:
// data
    size_t n_states;
    vector<SplitDistro> emis_distros;
    vector<double> start_trans;
    vector< vector<double> > trans;
    vector<double> end_trans;
};

HMM::HMM(const size_t n)
{
    set_states_n(n);
}

void
HMM::simulate(const size_t state,  pair<double, double>  &obs) const
{
    const string name_arg = emis_distros[state].tostring();
    const vector<string> name_split = smithlab::split_whitespace_quoted(name_arg);
    
    string name_arg_1, name_arg_2;
    if (name_split.size() == 3)
    {
        const string name = "pois";
        name_arg_1 = name + "," + name_split[1];
        name_arg_2 = name + " " + name_split[2];
    }
    else
    {
        const string name = "nbd";
        name_arg_1 = name + "," + name_split[1] + "," + name_split[2];
        name_arg_2 = name + "," + name_split[3] + "," + name_split[4];
    }

    Distro first_distro(name_arg_1), second_distro(name_arg_2);
    obs = std::make_pair(first_distro.sample(), second_distro.sample());
}

void 
HMM::simulate(const vector<size_t> &states,
              vector< pair<double, double> > &obs) const
{
    // get underlying distribution
    vector< pair<Distro, Distro> > distros;
    for (size_t i = 0; i < n_states; ++i)
    {
        const string name_arg = emis_distros[i].tostring();
        const vector<string> name_split = smithlab::split_whitespace_quoted(name_arg);
        
        string name_arg_1, name_arg_2;
        if (name_split.size() == 3)
        {
            const string name = "pois";
            name_arg_1 = name + "," + name_split[1];
            name_arg_2 = name + " " + name_split[2];
        }
        else
        {
            const string name = "nbd";
            name_arg_1 = name + "," + name_split[1] + "," + name_split[2];
            name_arg_2 = name + "," + name_split[3] + "," + name_split[4];
        }
        
        distros.push_back(std::make_pair(Distro(name_arg_1), Distro(name_arg_2)));
    }
    
    obs.clear();
    for (size_t i = 0; i < states.size(); ++i)
    {
        const size_t s = states[i];
        const pair<double, double> dp
            = std::make_pair(distros[s].first.sample(),
                             distros[s].second.sample());
        obs.push_back(dp);
    }
}

size_t 
get_state_label(const vector<double> &probs,
                Runif &runif)
{
    double p = runif.runif(0.0, 1.0);

    size_t i = 0;
    double s = probs.front();
    while (p > s) 
    {
        ++i;
        s += probs[i];
    }
    
    return i;
}

void 
HMM::simulate(const size_t n, vector<size_t> &states,
              vector< pair<double, double> > &obs) const
{   
    states.clear();
    Runif runif;
    states.push_back(get_state_label(start_trans, runif));
    for (size_t i = 1; i < n; ++i)
        states.push_back(get_state_label(trans[states[i-1]], runif));
    simulate(states, obs);
}

double
HMM::get_log_likelihood(const size_t &state,
                        const  pair<double, double>  &obs) const
{
    assert(state >= 0 && state < n_states);
    return emis_distros[state].log_likelihood(obs.first - obs.second);
}

double
HMM::get_log_likelihood(const size_t &state,
                        const vector< pair<double, double> > &obs) const
{
    assert(state >= 0 && state < n_states);
    vector<double> vals(obs.size());
    
    return emis_distros[state].log_likelihood(vals);
}

double
HMM::get_log_likelihood(const vector<size_t> &states,
                        const vector< pair<double, double> > &obs) const
{
    double llh = 0;
    assert(states.size() == obs.size());
    vector<double> vals(obs.size());
    
    for (size_t i = 0; i < states.size(); ++i)
    {
        assert(states[i] >= 0 && states[i] < n_states);
        vals[i] = obs[i].first - obs[i].second;
    }
    
    
    vector<double> log_start_trans(start_trans.size());
    for (size_t i = 0; i < start_trans.size(); ++i)
        log_start_trans[i] = log(start_trans[i]);

    vector< vector<double> > log_trans(trans.size(),
                                       vector<double>(trans.size(), 0));
    for (size_t i = 0; i < trans.size(); ++i)
        for (size_t j = 0; j < trans[i].size(); ++j)
            log_trans[i][j] = log(trans[i][j]);

    vector<double> log_end_trans(end_trans.size());
    for (size_t i = 0; i < end_trans.size(); ++i)
        log_end_trans[i] = log(end_trans[i]);
    
    llh += log_start_trans[states.front()] +
        emis_distros[states.front()].log_likelihood(vals.front());
    
    for (size_t i = 1; i < states.size(); ++i)
        llh += log_trans[states[i-1]][states[i]] +
            emis_distros[states[i]].log_likelihood(vals[i]);
    
     llh += log_end_trans[states.back()];
     
     return llh;
}

void
HMM::set_states_n(const size_t n)
{
    n_states = n;
    emis_distros.resize(n_states, SplitDistro("skel"));
    start_trans.resize(n_states);
    trans.resize(n_states);
    for (size_t i = 0; i < n_states; ++i)
        trans[i].resize(n_states);
    end_trans.resize(n_states);
}

static void
skip_comment_and_empty_lines(std::istream &stream, std::stringstream &ss)
{
    while (!stream.eof())
    {
        string str;
        getline(stream, str);
        const size_t comment_start = str.find("#");
        if (comment_start != string::npos)
            str.erase(comment_start);
        str = smithlab::strip(str);
        if (!str.empty())
            ss << str << endl;
    }
}

std::istream &operator>>(std::istream &in, HMM &hmm)
{
    std::stringstream ss(std::stringstream::in | std::stringstream::out);
    skip_comment_and_empty_lines(in, ss);
    
    size_t n(0);
    ss >> n; 
    string tmp_str;
    getline(ss, tmp_str);


    hmm.set_states_n(n);
    vector<SplitDistro> distros(n, SplitDistro("skel"));
    for (size_t i = 0; i < n; ++i)
    {
        string tmp_str;
        getline(ss, tmp_str);
        distros[i] = SplitDistro(tmp_str);
    }
    hmm.set_emission_distros(distros);
    
    vector<double> start_trans(n);
    for (size_t i = 0; i < n; ++i)
        ss >> start_trans[i];
    hmm.set_start_trans(start_trans);
    
    vector< vector<double> > trans(n, vector<double>(n, 0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            ss >> trans[i][j];
    hmm.set_trans(trans);

    vector<double> end_trans(n);
    for (size_t i = 0; i < n; ++i)
        ss >> end_trans[i];
    hmm.set_end_trans(end_trans);
    
    return in;
}

std::ostream &operator<<(std::ostream &out, const HMM &hmm)
{
    const size_t n = hmm.get_states_n();
    out << "# number of states" << endl;
    out << n << endl;
    
    out << "# emmission distributions" << endl;
    vector<SplitDistro> distros(hmm.get_emission_distros());
    copy(distros.begin(), distros.end(),
         ostream_iterator<SplitDistro>(out, "\n"));
    
    out << "# start to state transition probabilities" << endl;
    vector<double> start_trans(hmm.get_start_trans());
    copy(start_trans.begin(), start_trans.end(),
         ostream_iterator<double>(out, "\n"));

    out << "# state transition probabilities" << endl;
    vector< vector<double> > trans(hmm.get_trans());
    for (size_t i = 0; i < n; ++i)
    {
        copy(trans[i].begin(), trans[i].end(),
             ostream_iterator<double>(out, "\t"));
        out << endl;
    }
    
    out << "# states to end transition probabilities" << endl;
    vector<double> end_trans(hmm.get_end_trans());
    copy(end_trans.begin(), end_trans.end(),
         ostream_iterator<double>(out, "\n"));
    
    return out;
}

#endif
