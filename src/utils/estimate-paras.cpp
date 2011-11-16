/* estimate_paras.cpp
 * Song Qiang <qiang.song@usc.edu> 2010
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>


#include <cmath>
#include <cassert>

#include "OptionParser.hpp"
#include "Distro.hpp"
#include "SplitDistro.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using namespace std;


template <class T> void
read_data_file(const string &file_name,
               vector<T> &vals_a,
               vector<T> &vals_b)
{
    std::ifstream inf(file_name.c_str());
    T a, b;
    
    while ((inf >> a >> b) && !inf.bad())
    {
        vals_a.push_back(a);
        vals_b.push_back(b);
    }
}


template <class T> void
read_data_file(const string &file_name,
               vector<T> &vals)
{
    std::ifstream inf(file_name.c_str());
    T a;
    
    while((inf >> a) && !inf.bad())
    {
        vals.push_back(a);
    }
}

template <class Distro_Class> void
chi_square_test(const vector<double> &vals,
                const Distro_Class &distro,
                double &chi_square,
                size_t &df)
{
    const double tolerance = 1e-10;
    const size_t n = vals.size();

    // chi_square test requires that for each cell
    // there are at least 5 data points
    const double min_prob = 5 / static_cast<double>(n);
    
    // get the lower limit of bins
    double lower = 0;
    double current_prob = exp(distro.log_likelihood(lower));
    double next_prob = exp(distro.log_likelihood(lower - 1));
    while ((next_prob >= current_prob)
           || (next_prob < current_prob && next_prob > min_prob) )
    {
        lower -= 1;
        current_prob = next_prob;
        next_prob = exp(distro.log_likelihood(lower - 1));
    }
    
    double low_prob = 0;
    double tmp_val = lower - 1;
    double tmp_prob = exp(distro.log_likelihood(tmp_val));
    while (tmp_prob >= tolerance)
    {
        low_prob += tmp_prob;
        tmp_val -= 1;
        tmp_prob = exp(distro.log_likelihood(tmp_val));
    }
    
    // get the upper limit of bins
    double upper = 0;
    current_prob = exp(distro.log_likelihood(upper));
    next_prob = exp(distro.log_likelihood(upper + 1));
    while ((next_prob >= current_prob)
           || (next_prob < current_prob && next_prob > min_prob) )
    {
        upper += 1;
        current_prob = next_prob;
        next_prob = exp(distro.log_likelihood(upper + 1));
    }

    double up_prob = 0;
    tmp_val  = upper + 1;
    tmp_prob = exp(distro.log_likelihood(tmp_val));
    while (tmp_prob >= tolerance)
    {
        up_prob += tmp_prob;
        tmp_val += 1;
        tmp_prob = exp(distro.log_likelihood(tmp_val));
    }


    // set bins and count observed count
    vector<size_t> counts(static_cast<size_t>(upper - lower + 1), 0);
    size_t low_count(0), up_count(0);
    const int offset = static_cast<int>(-lower);
    
    for (size_t i = 0; i < vals.size(); ++i)
        if (vals[i] < lower)
            ++low_count;
        else if (vals[i] > upper)
            ++up_count;
        else
            ++counts[static_cast<size_t>(vals[i] + offset)];

    // compute chi square statistic
    chi_square = 0;
    if (low_prob >= tolerance)
        chi_square += pow(low_count - low_prob * n, 2) / (low_prob * n);
    if (up_prob >= tolerance)
        chi_square += pow(up_count - up_prob * n, 2) / (up_prob * n);
    
    for (size_t i = 0; i < counts.size(); ++i)
    {
        const double v = static_cast<double>(i) - offset;
        const double prob = exp(distro.log_likelihood(v));
        const double expected = n * prob;
        chi_square += pow(counts[i] - expected, 2) / expected;
    }
    
    // determine degree of freedom
    size_t non_empty_cells = 0;
    if (low_count != 0) ++non_empty_cells;
    if (up_count != 0) ++non_empty_cells;
    for (size_t i = 0; i < counts.size(); ++i)
        if (counts[i] != 0) ++non_empty_cells;

    const string str = distro.tostring();
    if (str.find("skel") != string::npos)
        df = non_empty_cells - (2 + 1);
    else if (str.find("nbdiff") != string::npos)
        df = non_empty_cells - (4 + 1);
    else if (str.find("pois") != string::npos)
        df = non_empty_cells - (1 + 1);
    else if (str.find("nbd") != string::npos)
        df = non_empty_cells - (2 + 1);
}

double
op_substract(double x, double y) 
{
    return x - y;
}

int
main(int argc, const char **argv)
{
    string distro_name;
    bool Goodness_of_Fit_Test = false;
    bool Log_Likelihood_Request = false;
    bool PMF_Table_Request = false;
    bool Model_Select_Criterion = false;

    OptionParser opt_parse("estimate_paras",
                           "This program estimate the parameters given data file",
                           "data_file");
    opt_parse.add_opt("distro", 'd', "Name of fitted distribution",
                      true, distro_name);
    opt_parse.add_opt("Goodness_of_Fit_Test", 't',
                      "Run Chi-Square goodness of fit test",
                      false, Goodness_of_Fit_Test);
    opt_parse.add_opt("Log_Likelihood_Request", 'l',
                      "Compute log likelihood",
                      false, Log_Likelihood_Request);
    opt_parse.add_opt("PMF_Table_Request", 'p',
                      "Output table of probabilties",
                      false, PMF_Table_Request);
    opt_parse.add_opt("Model_Select_Criterion", 'm',
                      "Compute AIC and BIC",
                      false, Model_Select_Criterion);

    vector<string> leftover_args;
    try
    {
        opt_parse.parse(argc, argv, leftover_args);
    }
    catch(RMAPOptionException e)
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

    try
    {
        
    const string file_name = leftover_args.front();

    if (distro_name == "skel" || distro_name == "nbdiff"
        || distro_name == "gauss" )
    {
        vector<double> vals_a, vals_b;
        read_data_file(file_name, vals_a, vals_b);

        assert(vals_a.size() == vals_b.size());
    
        vector<double> vals(vals_a.size());
        transform(vals_a.begin(), vals_a.end(),
                  vals_b.begin(), vals.begin(),
                  op_substract);

        const double max_val = *max_element(vals.begin(), vals.end());
        const double min_val = *min_element(vals.begin(), vals.end());

        SplitDistro distro(distro_name);
        distro.estimate_params_ml(vals_a, vals_b);

        cout << distro << "\t";
    
        if (Log_Likelihood_Request)
        {
            double llr = distro.log_likelihood(vals);
            cout << "Log Likelihood = " << llr << "\t";
        }

        if (Model_Select_Criterion)
        {
            double llr(0);
            for (size_t i = 0; i < vals.size(); ++i)
                llr += distro.log_likelihood(vals[i]);

            size_t k = 0;
            const string str = distro.tostring();
            if (str.find("skel") != string::npos)
                k = 2;
            else if (str.find("nbdiff") != string::npos)
                k = 4;
            else if (str.find("pois") != string::npos)
                k = 1;
            else if (str.find("nbd") != string::npos)
                k = 2;
            else if (str.find("gauss") != string::npos)
                k = 2;

            
            double AIC = 2 * k - 2 * llr;
            double BIC = log(vals.size()) * k - 2 * llr;
            cout << "AIC = " << AIC
                 << ", BIC = " << BIC << "\t";
        }
    
        if (Goodness_of_Fit_Test)
        {
            double chi_square(0);
            size_t df(0);
        
            chi_square_test(vals, distro, chi_square, df);
            cout << "Chi-Squared Test Statistic = " << chi_square
                 << ", DF = " << df << "\t";
        }
    
        cout << endl;

        if (PMF_Table_Request)
        {
            for (double r = min_val; r < max_val; ++r)
                cout << r << "\t"
                     << exp(distro.log_likelihood(r)) << endl;
        }
    }
    else if (distro_name == "pois" || distro_name == "nbd")
    {
        vector<double> vals;
        read_data_file(file_name, vals);

        const double max_val = *max_element(vals.begin(), vals.end());
        const double min_val = *min_element(vals.begin(), vals.end());

        Distro distro(distro_name);
        distro.estimate_params_ml(vals);

        cout << distro << "\t";
    
        if (Log_Likelihood_Request)
        {
            double llr = distro.log_likelihood(vals);
            cout << "Log Likelihood = " << llr << "\t";
        }

        if (Model_Select_Criterion)
        {
            double llr(0);
            for (size_t i = 0; i < vals.size(); ++i)
                llr += distro.log_likelihood(vals[i]);

            size_t k = 0;
            const string str = distro.tostring();
            if (str.find("skel") != string::npos)
                k = 2;
            else if (str.find("nbdiff") != string::npos)
                k = 4;
            else if (str.find("pois") != string::npos)
                k = 1;
            else if (str.find("nbd") != string::npos)
                k = 2;
            
            double AIC = 2*k - 2 * llr;
            double BIC = log(vals.size()) * k - 2 * llr;
            cout << "AIC = " << AIC
                 << ", BIC = " << BIC << "\t";

        }
    
        if (Goodness_of_Fit_Test)
        {
            double chi_square(0);
            size_t df(0);
        
            chi_square_test(vals, distro, chi_square, df);
            cout << "Chi-Squared Test Statistic = " << chi_square
                 << ", DF = " << df << "\t";
        }
    
        cout << endl;

        if (PMF_Table_Request)
        {
            for (double r = min_val; r < max_val; ++r)
                cout << r << "\t"
                     << exp(distro.log_likelihood(r)) << endl;
        }

    }
    else
    {
        cout << "unrecognized distribution name" << endl;
    }
    }
    catch (SMITHLABException e)
    {
        cerr << e.what() << endl;
    }

    return 0;
}
