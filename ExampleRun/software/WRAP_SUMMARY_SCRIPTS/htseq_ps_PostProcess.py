#! /bin/python3.9.5

import sys
import numpy as np
import argparse
import os
import pandas as pd
from functools import reduce
import re
import subprocess
from itertools import chain


####----Take and Parse Arguments----####

parser=argparse.ArgumentParser(description='Combine htseq, splicing deficiency, and PSI/PSO files to generate sample summary.')

parser.add_argument('-ht','--htseq', const=1, nargs='?', required=False, default=1, metavar='',help='path to htseq file')
parser.add_argument('-p','--psipso', const=1, nargs='?', required=False, default=1, metavar='', help='path to psi/pso file')
parser.add_argument('-oh', '--outht', type=argparse.FileType('w'), nargs='?', required=False, default='htseq_psProcessed.txt', metavar='', help='Processed HTSEQ outfile name')
parser.add_argument('-op', '--outps', type=argparse.FileType('w'), nargs='?', required=False, default='PS_htseqProcessed.txt', metavar='', help='Processed PSI/PSO outfile name')
parser.add_argument('-om', '--outmg', type=argparse.FileType('w'), nargs='?', required=False, default='merged_PS_htseqProcessed.txt', metavar='', help='Processed PSI/PSO with HTSEQ merged outfile name')


args=parser.parse_args()


####----Assign File Input----####

htseq=args.htseq
psipso=args.psipso
outh=args.outht
outp=args.outps
outm=args.outmg


ht = pd.read_csv(htseq, sep = '\t', header = 0)
ps = pd.read_csv(psipso, sep = '\t', header = 0)

## Remove Filter column if it exists
if ps.columns[-1] == "Filter":
	ps = ps.iloc[:,:-1]

## Split Exon ID to extract gene symbol to column
ps[["Gene"]] = ps["ExonID"].str.split('.', n = 1, expand = True)[0]

## add suffix to represent htseq and psi/pso data
ps = ps.add_suffix("_ps")
ps.columns = [ps.columns[:-1],'Gene']
ps.rename(columns = { ps.columns[0]: "ExonID"}, inplace = True)
ht = ht.add_suffix("_ht")
ht.rename(columns = { ht.columns[0]: "Gene"}, inplace = True)

## Merge data base on Gene Symbol
merge = pd.merge(ps, ht, on = ["Gene"], how = "inner")

## Generate lists of samples to split by
ht_cols = ['Gene','ExonID']
ht_cols.append(list(merge.filter(regex = '_ht$').columns))
ps_cols = ['Gene','ExonID',list(merge.filter(regex = '_ps$').columns)]
ht_samp = list(merge.filter(regex = '_ht$').columns)
ps_samp = list(merge.filter(regex = '_ps$').columns)

#separate sample types
ht_sub = merge.filter(regex = 'Gene|ExonID|_ht$')
ps_sub = merge.filter(regex = 'Gene|ExonID|_ps$')

## Move gene column to the first position
#first_col_mg = merge.pop("Gene")
#first_col_ht = ht_sub.pop("Gene")
#first_col_ps = ps_sub.pop("Gene")
#merge.insert(0, "Gene", first_col_mg)
#ht_sub.insert(0, "Gene", first_col_ht)
#ps_sub.insert(0, "Gene", first_col_ps)

## Remove Gene Column
merge.pop("Gene")
ht_sub.pop("Gene")
ps_sub.pop("Gene")

#Remove suffix
ht_sub.columns = ht_sub.columns.str.rstrip('_ht')
ps_sub.columns = ps_sub.columns.str.rstrip('_ps')




## Write to outfile
merge.to_csv(outm, sep = '\t', index = False, na_rep = "NaN")
ht_sub.to_csv(outh, sep = '\t', index = False, na_rep = "NaN")
ps_sub.to_csv(outp, sep = '\t', index = False, na_rep = "NaN")





