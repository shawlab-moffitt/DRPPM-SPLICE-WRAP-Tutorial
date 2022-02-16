#! /bin/python3.9.5

import sys
import numpy as np
import argparse

####----Compatibility----####

#compatible with RSeQC v4.0.0

#must include 7 input files in correct order
#sample input: python summarygen.py {SAMPLENAME}.Aligned.sortedByCoord.out.summary.txt {SAMPLENAME}_junction_annotation_summary_more.txt rseqc_bam_stat_report.txt {SAMPLENAME}.STAR.Log.final.out {SAMPLENAME}_infer_experiment.txt {SAMPLENAME}_inner_distance.txt {SAMPLENAME}_read_distribution.txt {SAMPLENAME}_intron_summary.txt
#tin.py                 -- {SAMPLENAME}.Aligned.sortedByCoord.out.summary.txt
#junction_annotation.py -- {SAMPLENAME}_junction_annotation_summary_more.txt
#bam_stat.py            -- rseqc_bam_stat_report.txt (sample specific)
#STAR_log_final         -- {SAMPLENAME}.STAR.Log.final.out (from STAR mapping)
#infer_experiment.py    -- {SAMPLENAME}_infer_experiment.txt
#inner_distance.py      -- {SAMPLENAME}_inner_distance.txt
#read_distribution.py   -- {SAMPLENAME}_read_distribution.txt
#splicing deficiency    -- {SAMPLENAME}_intron_summary.txt



parser=argparse.ArgumentParser(description='Combine summary files to generate global sample summary.')

parser.add_argument('-t','--tin', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='',help='path to TIN summary file')
parser.add_argument('-j','--juncanno', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to junction annotation summary file')
parser.add_argument('-b','--bamstat', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to bam stat summary file')
parser.add_argument('-l','--logfin', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to STAR log final output file')
parser.add_argument('-e','--inferexp', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to infer experiment summary file')
parser.add_argument('-d','--innerdist', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to inner distance summary file')
parser.add_argument('-r','--readdist', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to read distribution summary file')
parser.add_argument('-n','--intron', type=argparse.FileType('r'), const=1, nargs='?', required=False, default=1, metavar='', help='path to intron summary file')
parser.add_argument('-s','--samplename', required=False, nargs='?', default=1, metavar='', help='Sample name input to identify output file')
parser.add_argument('-R','--row', action='store_true', required=False, help='Outfile given in row format, if not specified both row and column file formats given')
parser.add_argument('-C','--column', action='store_true', required=False, help='Outfile given in column format, if not specified both row and column file formats given')

args=parser.parse_args()

####----Input----####
tin=args.tin
juncanno=args.juncanno
bamstat=args.bamstat
logfin=args.logfin
inferexp=args.inferexp
innerdist=args.innerdist
readdist=args.readdist
intron=args.intron
samp=args.samplename


mainheader=np.array(['init'])
mainstats=np.array(['init'])

####----TIN file----####
if tin == 1:
	print("TIN summary file missing")
else:
	tinheader=tin.read().splitlines()[0] #extract header
	tinheader=tinheader.split('\t') #split by tab delim
	TINmean=tinheader[1]
	TINmed=tinheader[2]
	TINsd=tinheader[3]
	mainheader=np.append(mainheader, ['Sample',TINmean,TINmed,TINsd])
	tin.seek(0)
	tinstat=tin.read().splitlines()[1]
	tinstat2=tinstat.split('\t')
	sampname=tinstat2[0].split('.')[0]
	TINmeannum=tinstat2[1]
	TINmednum=tinstat2[2]
	TINsdnum=tinstat2[3]
	mainstats=np.append(mainstats, [sampname,TINmeannum,TINmednum,TINsdnum])



####----junction annotation----####
found_sect = False
if juncanno == 1:
	print("Junction Annotation summary file missing")
else:
	for line in juncanno.read().splitlines():
		if line.startswith('='):
			found_sect = True
		elif found_sect:
			if ':' in line:
				head=line.replace(' ','_').split('\t')[0].split(':')[0]
				val=line.split('\t')[-1]
				mainheader=np.append(mainheader,[head])
				mainstats=np.append(mainstats,[val])



####----BAM stat report----####
if bamstat == 1:
	print("BAM Stat Report summary file missing")
else:
	for line in bamstat.read().splitlines()[5:]:
		if ':' in line:
			head=line.replace(' ','_').split(':')[0]
			if '<' in head:
				head='mapq_non_unique_reads'
			if '>=' in head:
				head='mapq_unique_reads'
			val=line.replace(' ','').split(':')[-1]
			mainheader=np.append(mainheader,[head])
			mainstats=np.append(mainstats,[val])
		elif 'Non primary hits' in line: #special instance for missing semicolon in file
			head='_'.join(line.split(' ')[0:3])
			val=line.split(' ')[-1]
			mainheader=np.append(mainheader,[head])
			mainstats=np.append(mainstats,[val])

# adjust for special characters
mainheader=[s.replace("'+'","pos") for s in mainheader] #positive
mainheader=[s.replace("'-'","neg") for s in mainheader] #negative



####----STAR log final out report----####
if logfin == 1:
	print("STAR Log Final Report summary file missing")
else:
	for line in logfin.read().splitlines()[5:] :
		if '|' in line:
			head=line.lstrip().replace(' ','_').split('|')[0]
			if '%' in head:
				head=head.replace('%','_fraction_')
			val=line.replace('\t','').split('|')[-1]
			mainheader=np.append(mainheader,[head])
			mainstats=np.append(mainstats,[val])



####----infer experiment----####
if inferexp == 1:
	print("Infer Experiment summary file missing")
else:
	failed=inferexp.read().splitlines()[3]
	failh=failed.replace(' ','_').split(':')[0]
	failv=failed.replace(' ','').split(':')[-1]
	inferexp.seek(0)
	forw=inferexp.read().splitlines()[4]
	frowh='_'.join(forw.split(' ')[0:4])+'_stranded_forward'
	frowv=forw.replace(' ','').split(':')[-1]
	inferexp.seek(0)
	rev=inferexp.read().splitlines()[5]
	revh='_'.join(rev.split(' ')[0:4])+'_stranded_reverse'
	revv=rev.replace(' ','').split(':')[-1]
	mainheader=np.append(mainheader,[failh,frowh,revh])
	mainstats=np.append(mainstats,[failv,frowv,revv])



####----inner distance----####
if innerdist == 1:
	print("Inner Distance summary file missing")
else:
	indis=innerdist.read().splitlines()[1]
	indis=indis.split('\t')
	sampname2='_'.join(indis[0].split('_')[0:5])
	indismean=indis[1]
	indismed=indis[2]
	indissd=indis[3]
	mainheader=np.append(mainheader,["inner_distance_mean","inner_distance_median","inner_distance_sd"])
	mainstats=np.append(mainstats,[indismean,indismed,indissd])



####----read distribution----####
if readdist == 1:
	print("Read Distribution summary file missing")
else:
	for line in readdist.read().splitlines()[5:14]:
		headtot=line.split()[0]+'_total_bases'
		valtot=line.split()[1]
		headtag=line.split()[0]+'_tags/Kb'
		valtag=line.split()[3]
		mainheader=np.append(mainheader,[headtot,headtag])
		mainstats=np.append(mainstats,[valtot,valtag])



####----intron summary----####
#if hasattr(args, 'intron') == True:
if intron == 1:
	print("Intron summary file missing")
else:
	intrh=intron.read().splitlines()[0]
	intrh=intrh.replace(' ','_').split('\t')
	intrh1n=intrh[1].replace('%','Fraction_')
	mainheader=np.append(mainheader,[intrh1n,intrh[2],intrh[3],intrh[4],intrh[5]])
	intron.seek(0)
	intrv=intron.read().splitlines()[1]
	intrv=intrv.split('\t')
	intrvperc=float(intrv[1])/100
	mainstats=np.append(mainstats,[intrvperc,intrv[2],intrv[3],intrv[4],intrv[5]])
	sampname3=intrv[0].split('.')[0]



####----remove/replace certain special characters----####
mainheader=[s.replace("'","_") for s in mainheader] #safe
mainheader=[s.replace("(","_") for s in mainheader] #safe
mainheader=[s.replace(")","") for s in mainheader]  #safe
mainheader=[s.replace("/","_") for s in mainheader] #safe
mainheader=[s.replace(":","_") for s in mainheader] #safe
mainheader=[s.replace("\"","_") for s in mainheader] #safe
mainheader=[s.replace("-","_") for s in mainheader] #safe
mainheader=[s.replace(",","_") for s in mainheader] #safe

## checks for no values and zeros and percentages turned to fractions
for i,j in enumerate(mainstats):
	if len(j) == 0:
		mainstats[i]='NA'
	if j == '0':
		mainstats[i]='NA'
	if '%' in j:
		x=float(mainstats[i].replace('%',''))
		mainstats[i]=x/100

## checks for missing items and zeros, removes leading and trailing '_' and capitalize first letter
for i,j in enumerate(mainheader):
	if '_' in j:
		mainheader[i]=mainheader[i].strip('_')
	mainheader[i]=mainheader[i].capitalize()
	if len(j) == 0:
		mainheader[i]='NA'
	if j == '0':
		mainheader[i]='NA'

mainheader=np.delete(mainheader, [0])
mainstats=np.delete(mainstats, [0])


####----Output----####
if samp == 1: #user did not give sample name
	if 'sampname' in locals(): #get samplename from TIN file
		if args.row:
			outfile=open("".join(sampname+"_summary_row.tsv"), 'w')
		if args.column:
			outfile2=open("".join(sampname+"_summary_col.tsv"), 'w')
		if (args.row == False) and (args.column == False):
			outfile=open("".join(sampname+"_summary_row.tsv"), 'w')
			outfile2=open("".join(sampname+"_summary_col.tsv"), 'w')
	elif 'sampname2' in locals(): #get samplename from inner distance file
		if args.row:
			outfile=open("".join(sampname2+"_summary_row.tsv"), 'w')
		if args.column:
			outfile2=open("".join(sampname2+"_summary_col.tsv"), 'w')
		if (args.row == False) and (args.column == False):
			outfile=open("".join(sampname2+"_summary_row.tsv"), 'w')
			outfile2=open("".join(sampname2+"_summary_col.tsv"), 'w')
	elif 'sampname3' in locals(): #get samplename from intron file
		if args.row:
			outfile=open("".join(sampname3+"_summary_row.tsv"), 'w')
		if args.column:
			outfile2=open("".join(sampname3+"_summary_col.tsv"), 'w')
		if (args.row == False) and (args.column == False):
			outfile=open("".join(sampname3+"_summary_row.tsv"), 'w')
			outfile2=open("".join(sampname3+"_summary_col.tsv"), 'w')
	else: #generic sample name
		if args.row:
			outfile=open("sample_summary_row.tsv", 'w')
		if args.column:
			outfile2=open("sample_summary_col.tsv", 'w')
		if (args.row == False) and (args.column == False):
			outfile=open("sample_summary_row.tsv", 'w')
			outfile2=open("sample_summary_col.tsv", 'w')
else: #user given sample name
	if args.row:
		outfile=open("".join(samp+"_summary_row.tsv"), 'w')
	if args.column:
		outfile2=open("".join(samp+"_summary_col.tsv"), 'w')
	if (args.row == False) and (args.column == False):
		outfile=open("".join(samp+"_summary_row.tsv"), 'w')
		outfile2=open("".join(samp+"_summary_col.tsv"), 'w')

## Write numpy array to tab delimited file
#as rows
if args.row:
	summarray=np.vstack((mainheader,mainstats))
	np.savetxt(outfile, summarray, fmt='%s', delimiter='\t')
#as columns
if args.column:
	summarray2=np.column_stack((mainheader,mainstats))
	np.savetxt(outfile2, summarray2, fmt='%s', delimiter='\t')
#both
if (args.row == False) and (args.column == False):
	summarray=np.vstack((mainheader,mainstats))
	np.savetxt(outfile, summarray, fmt='%s', delimiter='\t')
	summarray2=np.column_stack((mainheader,mainstats))
	np.savetxt(outfile2, summarray2, fmt='%s', delimiter='\t')

