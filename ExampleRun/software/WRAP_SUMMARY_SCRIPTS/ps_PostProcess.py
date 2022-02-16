#! /bin/python3.9.5

import sys
import numpy as np
import argparse
import os
import pandas as pd
from functools import reduce
import re
import subprocess


####----Take and Parse Arguments----####

parser=argparse.ArgumentParser(description='Filter PSI/PSO files.')

parser.add_argument('-p','--psipso', const=1, nargs='?', required=False, default=1, metavar='', help='path to psi/pso file')
parser.add_argument('-o', '--output', required=True, metavar='', help='Outfile name with path if desired (No extension)')


args=parser.parse_args()

####----Assign File Input----####

input_file=args.psipso
outfile=args.output

out_snf = open(''.join(outfile+"_Summary_unfiltered.txt"), 'w')
out_sf = open(''.join(outfile+"_Summary_filtered.txt"), 'w')
out_f = open(''.join(outfile+"_filter.txt"), 'w')


with open(input_file) as ps:
	for line in ps:

		##--Start Counter--##

		one = 0
		zero = 0
		NA = 0

		##--Read line by line--##

		line = line.strip('\n').split('\t')

		##--Write header to files--##

		if line[0] == "ExonID":
			line.append("Filter")
			out_snf.write('\t'.join(line)+'\n')
			out_sf.write('\t'.join(line)+'\n')
			fil = [line[0],"Filter"]
			out_f.write('\t'.join(fil)+'\n')

		##--Evaluate contents of file--##

		else:
			#Get number of samples/columns
			sampnum = len(line) - 1

			##--Count number of 1, 0, or NAs in line--##

			for i in line:
				if i == '1.0':
					one = one + 1
				elif i == '0.0':
					zero = zero + 1
				elif i == 'NA':
					NA = NA + 1
				else:
					continue

			##--Write out lines to respective files--##

			if one + NA == sampnum:
				line.append("Remove")
				out_snf.write('\t'.join(line)+'\n')
				fil = [line[0],"Remove"]
				out_f.write('\t'.join(fil)+'\n')
			elif zero + NA == sampnum:
				line.append("Remove")
				out_snf.write('\t'.join(line)+'\n')
				fil = [line[0],"Remove"]
				out_f.write('\t'.join(fil)+'\n')
			else:
				line.append("Keep")
				out_snf.write('\t'.join(line)+'\n')
				out_sf.write('\t'.join(line)+'\n')
				fil = [line[0],"Keep"]
				out_f.write('\t'.join(fil)+'\n')


out_snf.close()
out_sf.close()
out_f.close()


