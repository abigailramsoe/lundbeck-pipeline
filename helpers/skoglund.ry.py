#!/usr/bin/env python

"""
Version: 	0.4
Author:		Pontus Skoglund
Contact: 	pontus.skoglund@gmail.com
Date: 		28 June 2013
Citation: 	P. Skoglund, J.Stora, A. Gotherstrom, M. Jakobsson (2013)
		Accurate sex identification of ancient human remains using DNA shotgun sequencing.
		Journal of Archaeological Science

Usage:		python ry_compute.py <SAM formatted data from stdin>

Example:	samtools view -q 30 mybamfile.bam | python XYkaryotyper.py

		(for specification on the SAM format and a the samtools suite, see Li, Handsaker et al. 2009, Bioinformatics)

Output:		[Total number of alignments in input] [Number of X and Y alignments identified] [R_y] [R_y standard error] [95% CI for R_Y] [Inferred sex]
"""

import sys
import math
from optparse import OptionParser

usage = "usage: %prog [options] <SAM formatted data from stdin>"
parser = OptionParser(usage=usage)
parser.add_option("--chrXname", action="store", type="string", dest="chrXname",help="Identifier for the X chromosome in the SAM input (use if different than chrX, X etc)",default="X")
parser.add_option("--chrYname", action="store", type="string", dest="chrYname",help="Identifier for the Y chromosome in the SAM input (use if different than chrY, Y etc)",default="Y")
parser.add_option("--malelimit", action="store", type="float", dest="malelimit",help="Upper R_y limit for assignment as XY/male",default=0.075)
parser.add_option("--femalelimit", action="store", type="float", dest="femalelimit",help="Lower R_y limit for assignment as XX/female",default=0.016)
parser.add_option("--digits", action="store", type="int", dest="digits",help="Number of decimal digits in R_y output",default=4)
parser.add_option("--noheader", action="store_true", dest="noheader",help="Do not print header describing the columns in the output",default=False)
(options, args) = parser.parse_args()

chrYcount=0
chrXcount=0
totalcount=0
for line in sys.stdin:

	#SAM header lines are skipped
	if line[0] == '@':continue
	totalcount += 1

	col=line.split()
	chromosome=col[2]

	if options.chrYname in chromosome: chrYcount += 1
	elif options.chrXname in chromosome: chrXcount += 1


#compute R_y with 95% confidence interval
n=chrYcount+chrXcount # total number of used alignments

if n == 0:
	if options.noheader == False:
		print('Nseqs\tNchrY+NchrX\tNchrY\tR_y\tSE\t95% CI\tAssignment')
	print(totalcount,'\t',n,'\t',0,'\t',0,'\t',0,'\t','0-0','\t',"N/A")
	exit(0)

Ry=1.0*chrYcount/n
SE=math.sqrt((Ry*(1.0-Ry))/n)
confinterval=1.96*SE


#use criteria to infer chromosomal sex
gender='NA'
if (Ry < options.femalelimit) and (Ry > options.malelimit):
	gender='Not Assigned'
elif Ry==0.0:
	gender='consistent with XX'
elif Ry+confinterval < options.femalelimit:
	gender='XX'
elif Ry-confinterval > options.malelimit:
        gender='XY'
elif Ry-confinterval > options.femalelimit and Ry+confinterval > options.malelimit:
	gender='consistent with XY but not XX'
elif Ry-confinterval < options.femalelimit and Ry+confinterval < options.malelimit:
	gender='consistent with XX but not XY'
else:
	gender='Not Assigned'

if options.noheader == False:
	print('Nseqs\tNchrY+NchrX\tNchrY\tR_y\tSE\t95% CI\tAssignment')
print(totalcount,'\t',n,'\t',chrYcount,'\t',round(Ry,options.digits),'\t',round(SE,options.digits),'\t',str(round(Ry-confinterval,options.digits))+'-'+str(round(Ry+confinterval,options.digits)),'\t',gender)
