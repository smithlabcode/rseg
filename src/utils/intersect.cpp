/* intersect.cpp
 * Song Qiang <qiang.song@usc.edu> 2010
 */

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

bool
parse_bedline(const string &line,
			  string &chrom,
			  size_t &start,
			  size_t &end)
{
    stringstream ss(line);
    ss >> chrom >> start >> end;

    return !ss.fail();
}


int
main(int argc, char ** argv)
{
    if (argc < 3)
    {
        cerr << "intersect bedfile_1 bedfile_2" << endl;
        return 0;
    }

    ifstream infile_1(argv[1], fstream::in);	
    ifstream infile_2(argv[2], fstream::in);

    string chrom1, chrom2, line1, line2;
    size_t  start(0), end(0),
        start1(0), end1(0),
        start2(0), end2(0);	

    getline(infile_1, line1);
    bool flag1 = parse_bedline(line1, chrom1, start1, end1);
    getline(infile_2, line2);
    bool flag2 = parse_bedline(line2, chrom2, start2, end2);

    while (flag1 && flag2)
    {
        start = max(start1, start2);
        end = min(end1, end2);
        
        if (chrom1 == chrom2 && start < end)
        {
            cout << chrom1 << "\t"
                 << start << "\t"
                 << end << "\t"
                 << line1 << "\t"
                 << line2 << endl;
        }

        if (( chrom1 < chrom2 ) || (chrom1 == chrom2 && end1 < end2))
        {
            line1.clear();
            getline(infile_1, line1);
            flag1 = parse_bedline(line1, chrom1, start1, end1);
        } 
        else
        {
            line2.clear();
            getline(infile_2, line2);
            flag2 = parse_bedline(line2, chrom2, start2, end2);
        }
    }
    
    return 0;
}
