#! /bin/python3.9.5

import sys
import argparse
import re

parser=argparse.ArgumentParser(description='Get length annotation from GTF file.')

parser.add_argument('-g','--gtf', const=1, nargs='?', required=False, default=1, metavar='',help='Path to reference GTF file')
parser.add_argument('-o', '--output', required=True, metavar='', help='Outfile file name with path')

args=parser.parse_args()

gtf_file=args.gtf
outfile=args.output

out_file = open(outfile, 'w')

with open(gtf_file) as gtf:

	out_file.write('IntronID'+'\t'+'Length'+'\n')

	length_lst = []

	for line in gtf:
		line = line.strip('\n').split('\t')
		anno = line[8].split(';')
		anno_lst = [x.strip() for x in anno]

		#Gene ID
		gene_id = [anno_lst for anno_lst in anno_lst if 'gene_id ' in anno_lst]
		gene_id = gene_id.pop().split(' ')[1].strip('"')

		START = int(line[3])
		STOP = int(line[4])
		length = STOP - START

		pair = ''.join(gene_id+'\t'+str(length))

		length_lst.append(pair)

	uniq_pair_lst = list(set(length_lst))

	out_file.write('\n'.join(uniq_pair_lst))


out_file.close()