#! /usr/bin/python

# bed2html.py
# Song Qiang <qiang.song@usc.edu> 2010

import os, sys, math
from optparse import OptionParser

usage = "usage: bed2html [options] bed_file [html_file]"

parser = OptionParser()
parser.add_option("-p", "--proportion", type="float", dest="proportion", default = 0.618, 
				  help="Proportion of display covered by domain (default 0.618))")

(options, args) = parser.parse_args()

if len(args) == 0:
	parser.error("Need bed file name")
	
infile = args[0]
if len(args) == 1:
	outfile = os.path.splitext(infile)[0] + ".html"
else:
	outfile = args[1]

inf = open(infile)
outf = open(outfile, "w")
	
html_header = """
<html>
<head>
<title> </title>
<link rel="stylesheet" href="css/table.css" type="text/css">
</head>
<body>
Last updated: <!--#echo var="LAST_MODIFIED"\ -->
<br>
<table class="dmrtable" cellspacing="0">
"""
html_footer = """</table>
</body>
</html>
"""

table_odd_row_templ = """
<tr class="odd"> <td>%(ID)d</td><td><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%(CHROM)s%%3A%(START)d-%(END)d">%(CHROM)s:%(START)d-%(END)d</a></td>%(OTHER)s</tr>
"""

table_even_row_templ = """
<tr> <td>%(ID)d</td><td><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%(CHROM)s%%3A%(START)d-%(END)d">%(CHROM)s:%(START)d-%(END)d</a></td>%(OTHER)s</tr>
"""

outf.write(html_header + "\n")

id = 1

for line in inf:
	parts = line.split()
	if len(parts) < 3:
		break
	chrom = parts[0]
	start = int(parts[1])
	end = int(parts[2])
	width = end - start
	offset = math.floor(width * (1 / options.proportion - 1) / 2)
	start -= offset
	if start <= 0:
		start = 1
	end += offset
	other = ""
	for item in parts[3:]:
		other += "<td>" + item + "</td>"
	if id % 2 == 1:
		row = table_odd_row_templ % {"ID":id, \
									 "CHROM":chrom, \
									 "START":start, \
									 "END":end,\
									 "OTHER":other}
	else:
		row = table_even_row_templ % {"ID":id, \
									  "CHROM":chrom, \
									  "START":start, \
									  "END":end,\
									  "OTHER":other}
	outf.write(row)	
	
outf.write(html_footer + "\n")

inf.close()
outf.close()
