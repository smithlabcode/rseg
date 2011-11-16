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
    string chrom1, chrom, line1, line;
    size_t  start(0), end(0),
        start1(0), end1(0);

    getline(cin, line);
    bool flag = parse_bedline(line, chrom, start, end);

    while (flag)
    {
        line1.clear();
        getline(cin, line1);
        flag = parse_bedline(line1, chrom1, start1, end1);

        if (flag)
        {
            if ((chrom < chrom1)  ||
                (chrom == chrom1 && end < start1))
            {
                cout << chrom << "\t"
                     << start << "\t"
                     << end << "\t"
                     << line << endl;

                chrom = chrom1;
                start = start1;
                end = end1;
                line = line1;
            }
            else
            {
                end = max(end, end1);
                line += line1;
            }
        }
    }
    
    return 0;
}
