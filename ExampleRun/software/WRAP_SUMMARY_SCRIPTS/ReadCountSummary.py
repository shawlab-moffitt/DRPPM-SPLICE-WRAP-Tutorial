#! /bin/python3.6.8

import sys
import numpy as np
import argparse
import os
import pandas as pd
from functools import reduce
import re
import subprocess


####------------------------------------Script Function----------------------------------------####
# This script will take the input of the htseq_gene.lst, splicing_deficiency.lst, and psi_pso.lst #
# These files contain paths to files for all the samples                                          #
# It will take the information for all the samples files and merge that into a summary matrix     #
# Output will be a plain tab delim text file and HDF5 file                                        #
####-------------------------------------------------------------------------------------------####



####----Take and Parse Arguments----####

parser=argparse.ArgumentParser(description='Combine htseq, splicing deficiency, and PSI/PSO files to generate sample summary.')

parser.add_argument('-ht','--htseq', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='',help='path to htseq_gene_files.lst file')
parser.add_argument('-e','--exon', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='',help='path to htseq_exon_files.lst file')
parser.add_argument('-i','--intron', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='',help='path to htseq_intron_files.lst file')
parser.add_argument('-s','--splice', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to splicing_deficiency_files.lst file')
parser.add_argument('-p','--psipso', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to psi_pso_files.lst file')
#parser.add_argument('-j','--junct', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to juncsalvager.lst file')

args=parser.parse_args()


####----Assign File Input----####

htseq=args.htseq
exon=args.exon
intron=args.intron
splice=args.splice
psipso=args.psipso
#junct=args.junct


####----HTSEQ Gene File----####

if htseq == 1:
	print("HTSEQ file missing")
else:

	##--Blank lists to add dataframes to--##

	df_to_merge_counts = []
	df_to_merge_fpkm = []

	##--Read file line-by-line and extract needed file paths--##

	for line in htseq: #read file
		line = line.strip('\n').split('\t') #segment file lines into list of info
		sampName = line[0]    				#assign sample name
		fpkm_file = line[1]   				#assign fpkm file
		count_file = line[2]  				#assign rawcount file

		##--Read in count and FPKM files--##

		count = pd.read_csv(count_file, sep = '\t', header = 0, names = ['Gene',sampName]) 	  #read in rawcount file
		fpkm = pd.read_csv(fpkm_file, sep = '\t', header = 0, names = ['Gene',sampName])          #read in fpkm file

		##--Add to list to merge--##

		df_to_merge_counts.append(count)
		df_to_merge_fpkm.append(fpkm)

	##--Merge all dataframes in counts and fpkm list based on Gene column--##

	df_merged_counts = reduce(lambda  left,right: pd.merge(left,right,on=['Gene'],how='outer'), df_to_merge_counts).fillna('NA')
	df_merged_fpkm = reduce(lambda  left,right: pd.merge(left,right,on=['Gene'],how='outer'), df_to_merge_fpkm).fillna('NA')

	##--Remove decimal and number in ensemble ID--##

	df_merged_counts['Gene'] = df_merged_counts['Gene'].str.split('.').str.get(0)
	df_merged_fpkm['Gene'] = df_merged_fpkm['Gene'].str.split('.').str.get(0)


	##--Convert to Gene Symbol and Write Outfile--##

	#-Counts Summary File-#

	#Ensemble write outfile file
	df_merged_counts.to_csv(r'HTSEQ_Counts_Summary_Ensemble.txt', sep = '\t', index = False)
	#Convertion to gene symbol and write outfile
	bashCommand = 'drppm -CleanEnsemblGeneID2GeneName HTSEQ_Counts_Summary_Ensemble.txt /share/dept_bbsr/Projects/Shaw_Timothy/3352_Splicing_Pipeline_2021/references/rseqc/gencode.v22.annotation.gtf HTSEQ_Counts_Summary_Symbol.txt'
	process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
	output, error = process.communicate()


	#-FPKM Summary File-#

	#Ensemble write outfile file
	df_merged_fpkm.to_csv(r'HTSEQ_FPKM_Summary_Ensemble.txt', sep = '\t', index = False)
	#Convertion to gene symbol and write outfile
	bashCommand = 'drppm -CleanEnsemblGeneID2GeneName HTSEQ_FPKM_Summary_Ensemble.txt /share/dept_bbsr/Projects/Shaw_Timothy/3352_Splicing_Pipeline_2021/references/rseqc/gencode.v22.annotation.gtf HTSEQ_FPKM_Summary_Symbol.txt'
	process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
	output, error = process.communicate()


####----HTSEQ Exon Files----####

if exon == 1:
	print("HTSEQ Exon file missing")
else:
	##--Blank lists to add dataframes to--##

	df_to_merge_counts = []
	df_to_merge_fpkm = []

	##--Read file line-by-line and extract needed file paths--##

	for line in exon: #read file
		line = line.strip('\n').split('\t') #segment file lines into list of info
		sampName = line[0]    				#assign sample name
		fpkm_file = line[1]   				#assign fpkm file
		count_file = line[2]  				#assign rawcount file

		##--Read in count and FPKM files--##

		count = pd.read_csv(count_file, sep = '\t', header = 0, names = ['ExonID',sampName]) 	    #read in rawcount file
		fpkm = pd.read_csv(fpkm_file, sep = '\t', header = 0, names = ['ExonID',sampName])          #read in fpkm file

		##--Add to list to merge--##

		df_to_merge_counts.append(count)
		df_to_merge_fpkm.append(fpkm)

	##--Merge all dataframes in counts and fpkm list based on ExonID column--##

	df_merged_counts = reduce(lambda  left,right: pd.merge(left,right,on=['ExonID'],how='outer'), df_to_merge_counts).fillna('NA')
	df_merged_fpkm = reduce(lambda  left,right: pd.merge(left,right,on=['ExonID'],how='outer'), df_to_merge_fpkm).fillna('NA')

	##--Write Outfile--##

	#-Counts Summary File-#
	df_merged_counts.to_csv(r'HTSEQ_Exon_Counts_Summary.txt', sep = '\t', index = False)
	#-FPKM Summary File-#
	df_merged_fpkm.to_csv(r'HTSEQ_Exon_FPKM_Summary.txt', sep = '\t', index = False)


####----HTSEQ Intron Files----####

if exon == 1:
	print("HTSEQ Intron file missing")
else:
	##--Blank lists to add dataframes to--##

	df_to_merge_counts = []
	df_to_merge_fpkm = []

	##--Read file line-by-line and extract needed file paths--##

	for line in exon: #read file
		line = line.strip('\n').split('\t') #segment file lines into list of info
		sampName = line[0]    				#assign sample name
		fpkm_file = line[1]   				#assign fpkm file
		count_file = line[2]  				#assign rawcount file

		##--Read in count and FPKM files--##

		count = pd.read_csv(count_file, sep = '\t', header = 0, names = ['IntronID',sampName]) 	    #read in rawcount file
		fpkm = pd.read_csv(fpkm_file, sep = '\t', header = 0, names = ['IntronID',sampName])          #read in fpkm file

		##--Add to list to merge--##

		df_to_merge_counts.append(count)
		df_to_merge_fpkm.append(fpkm)

	##--Merge all dataframes in counts and fpkm list based on IntronID column--##

	df_merged_counts = reduce(lambda  left,right: pd.merge(left,right,on=['IntronID'],how='outer'), df_to_merge_counts).fillna('NA')
	df_merged_fpkm = reduce(lambda  left,right: pd.merge(left,right,on=['IntronID'],how='outer'), df_to_merge_fpkm).fillna('NA')

	##--Write Outfile--##

	#-Counts Summary File-#
	df_merged_counts.to_csv(r'HTSEQ_Intron_Counts_Summary.txt', sep = '\t', index = False)
	#-FPKM Summary File-#
	df_merged_fpkm.to_csv(r'HTSEQ_Intron_FPKM_Summary.txt', sep = '\t', index = False)


####----Splicing Deficiency File----####

if splice == 1:
	print("Splicing deficiency file missing")
else:
	df_to_merge_ens = [] #list to add SD dfs to
	df_to_merge_sym = [] #list to add SD dfs to

	##--Read file line-by-line and extract needed file paths--##

	for line in splice: #read file
		line = line.strip('\n').split('\t') #segment file lines into list of info
		sampName = line[0]    				#assign sample name
		SD_file_ens = line[1]
		SD_file_sym = line[2]  				    #assign SD with gene symbol file

		##--Read in SD files--##

		SD_ens = pd.read_csv(SD_file_ens, sep = '\t', header = 0, names = ['Gene',sampName], usecols = [0,1])     #read in fpkm file
		SD_sym = pd.read_csv(SD_file_sym, sep = '\t', header = 0, names = ['Gene',sampName], usecols = [0,1])     #read in fpkm file


		##--Add df to list--##

		df_to_merge_ens.append(SD_ens)
		df_to_merge_sym.append(SD_sym)

	##--Merge all dfs in list--##

	df_merged_SD_ens = reduce(lambda  left,right: pd.merge(left,right,on=['Gene'],how='outer'), df_to_merge_ens).fillna('NA')
	df_merged_SD_sym = reduce(lambda  left,right: pd.merge(left,right,on=['Gene'],how='outer'), df_to_merge_sym).fillna('NA')

	print(df_merged_SD_ens)
	print(df_merged_SD_sym)

	##--Write Splicing Deficiency Summary File--##

	df_merged_SD_ens.to_csv(r'Splicing_Deficiency_Summary_Ensemble.txt', sep = '\t', index = False)
	df_merged_SD_sym.to_csv(r'Splicing_Deficiency_Summary_Symbol.txt', sep = '\t', index = False)


####----PSI/PSO File----####

if psipso == 1:
	print("PSI/PSO file missing")
else:
	df_to_merge_psi = [] #list to add PSI dfs to
	df_to_merge_pso = [] #list to add PSO dfs to
	df_to_merge_3ASpsi = [] #list to add 3 prime AS_PSI dfs to
	df_to_merge_5ASpsi = [] #list to add 5 prime AS_PSI dfs to

	##--Read file line-by-line and extract needed file paths--##

	for line in psipso: #read file
		line = line.strip('\n').split('\t') #segment file lines into list of info
		sampName = line[0]    				#assign sample name
		psi_file = line[4]   				#assign PSI file
		pso_file = line[5]  				#assign PSO file
		psi_3p_file = line[6]				#assign 3 prime AltSplice PSI file
		psi_5p_file = line[7]				#assign 5 prime AltSplice PSI file


		##--Read in PSI/PSO files--##

		psi = pd.read_csv(psi_file, sep = '\t', header = 0, names = ['ExonID',sampName])                    #read in PSI file
		pso = pd.read_csv(pso_file, sep = '\t', header = 0, names = ['ExonID',sampName])                    #read in PSO file
		psi_3AS = pd.read_csv(psi_3p_file, sep = '\t', header = 0, names = ['ExonID',sampName])     #read in AltSplicPSI file
		psi_5AS = pd.read_csv(psi_5p_file, sep = '\t', header = 0, names = ['ExonID',sampName])     #read in AltSplicPSI file

		##--Add df to lists--##

		df_to_merge_psi.append(psi)
		df_to_merge_pso.append(pso)
		df_to_merge_3ASpsi.append(psi_3AS)
		df_to_merge_5ASpsi.append(psi_5AS)

	##--Merge based on exon ID--##

	df_merged_psi = reduce(lambda  left,right: pd.merge(left,right,on=['ExonID'],how='outer'), df_to_merge_psi).fillna('NA')
	df_merged_pso = reduce(lambda  left,right: pd.merge(left,right,on=['ExonID'],how='outer'), df_to_merge_pso).fillna('NA')
	df_merged_3ASpsi = reduce(lambda  left,right: pd.merge(left,right,on=['ExonID'],how='outer'), df_to_merge_3ASpsi).fillna('NA')
	df_merged_5ASpsi = reduce(lambda  left,right: pd.merge(left,right,on=['ExonID'],how='outer'), df_to_merge_5ASpsi).fillna('NA')


	##--Write PSI/PSO Summary Files--##

	#PSI
	df_merged_psi.to_csv(r'PSI_Summary.txt', sep = '\t', index = False, na_rep = "NaN")
	#PSO
	df_merged_pso.to_csv(r'PSO_Summary.txt', sep = '\t', index = False, na_rep = "NaN")
	#3 prime PSI
	df_merged_3ASpsi.to_csv(r'PSI_3_prime_alt_spice_Summary.txt', sep = '\t', index = False, na_rep = "NaN")
	#5 prime PSI
	df_merged_5ASpsi.to_csv(r'PSI_5_prime_alt_spice_Summary.txt', sep = '\t', index = False, na_rep = "NaN")



