#! /bin/python3.9.5

import sys
import argparse
import re

####----Take and Parse Arguments----####

parser=argparse.ArgumentParser(description='Convert BED file to GTF file.')

parser.add_argument('-b','--bed', const=1, nargs='?', required=False, default=1, metavar='',help='Path to BED file')
parser.add_argument('-g','--gtf', const=1, nargs='?', required=False, default=1, metavar='',help='Path to reference GTF file')
parser.add_argument('-s', '--source', required=True, metavar='', help='Source for column of GTF')
parser.add_argument('-f', '--feature', required=True, metavar='', help='Feature for column of GTF')
#parser.add_argument('-n', '--name', required=True, metavar='', help='Type of name from BED file for Attribute section in GTF')
parser.add_argument('-o', '--output', required=True, metavar='', help='Outfile file name with path')

args=parser.parse_args()


####----Assign File Input----####

bed_file = args.bed
gtf_file = args.gtf
sorc_b = args.source
feat_b = args.feature
#name_b = args.name
outfile = args.output

out_GTF = open(outfile, 'w')

####----Read Through GTF File----####

with open(gtf_file) as gtf:

	gene_dict = {}
	linenum=0

	for line in gtf:

		##--Read line by line--##

		if line.startswith("#"):
			continue
		else:
			line = line.strip('\n').split('\t')
			frame = line[7]
			anno = line[8].split(';')
			anno_lst = [x.strip() for x in anno]

			##--Assign variables--##

			#Gene ID
			gene_id = [anno_lst for anno_lst in anno_lst if 'gene_id ' in anno_lst]				#match with BED
			if len(gene_id) > 0:
				gene_id = gene_id.pop().split(' ')[1].strip('"')
				gene_id = gene_id.split('.')[0]
			else:
				continue
				#gene_id = 'NA'

			#Gene Type
			gene_tp = [anno_lst for anno_lst in anno_lst if 'gene_type ' in anno_lst]			#bring over
			if len(gene_tp) > 0:
				gene_tp = gene_tp.pop()
			else:
				continue
				#gene_tp = 'NA'

			#Gene Name
			gene_name = [anno_lst for anno_lst in anno_lst if 'gene_name ' in anno_lst]		#use to make new ID
			if len(gene_name) > 0:
				gene_name = gene_name.pop().split(' ')[1].strip('"')
			else:
				continue
				#gene_name = 'NA'

			#Transcript Type
			trans_type = [anno_lst for anno_lst in anno_lst if 'transcript_type ' in anno_lst]	#bring over
			if len(trans_type) > 0:
				trans_type = trans_type.pop()
			else:
				continue
				#trans_type = 'NA'

			#Transcript Name
			trans_name = [anno_lst for anno_lst in anno_lst if 'transcript_name ' in anno_lst]	#bring over
			if len(trans_name) > 0:
				trans_name = trans_name.pop()
			else:
				continue
				#trans_name = 'NA'

			#Exon Number
			int_num = [anno_lst for anno_lst in anno_lst if 'exon_number ' in anno_lst]			#bring over
			if len(int_num) > 0:
				#int_num = int_num.pop()
				int_num = int_num.pop().split(' ')[1]
				int_num = ''.join("intron_number "+int_num)
			else:
				continue
				#int_num = 'NA'

			#Level
			level = [anno_lst for anno_lst in anno_lst if re.search(r'level\s\d',anno_lst)]		#bring over
			if len(level) > 0:
				level = level.pop().split(' ')[1]
			else:
				continue
				#level = 'NA'

			##--Add to dictionary--##

			gene_dict[gene_id] = [gene_tp,gene_name,trans_type,trans_name,int_num,level,frame]


			linenum+=1
			if linenum%1000==0: #check status
				status=str(linenum)+' finished'
				sys.stderr.write(status+'\r')


####----Read Through BED File----####

with open(bed_file) as bed:

	linenum2 = 0

	for line in bed:

		##--Read line by line--##

		line = line.strip('\n').split('\t')

		##--Assign varaibles--##

		CHR = line[0]					#Chromosome
		START = line[1]					#Start Position
		STOP = line[2]					#Stop Position
		NAME = line[3].split('.')[0]	#Name
		SCORE = line[4]					#Score
		STR = line[5]					#Strand

		##--Get dictionary values based on Gene ID--##

		if NAME in gene_dict:
			info = gene_dict.get(NAME)
			gene_tp = info[0]
			gene_name = info[1]
			trans_type = info[2]
			trans_name = info[3]
			int_num = info[4]
			level = info[5]
			frame = info[6]

			##--Make new intron_id--##

			cont = [str(gene_name),str(CHR),str(START),str(STOP),str(STR)]
			ID = '_'.join(cont)
			#Apply to IDs
			g_id = "".join("gene_id \""+ID+"\"")
			g_name="".join("gene_name \""+ID+"\"")
			t_id = "".join("transcript_id \""+ID+"\"")
			i_id = "".join("intron_id \""+ID+"\"")

			##--Join old and new variables in correct order--##

			#int_num and/or level not always present, only write out if they exist
			if len(int_num) > 0 and len(level) > 0:
				attr_lst = [g_id,t_id,gene_tp,g_name,trans_type,trans_name,int_num,i_id,''.join("level "+str(level))]
			elif len(int_num) > 0 and len(level) == 0:
				attr_lst = [g_id,t_id,gene_tp,g_name,trans_type,trans_name,int_num,i_id]
			elif len(int_num) == 0 and len(level) > 0:
				attr_lst = [g_id,t_id,gene_tp,g_name,trans_type,trans_name,i_id,''.join("level "+str(level))]
			elif len(int_num) == 0 and len(level) == 0:
				attr_lst = [g_id,t_id,gene_tp,g_name,trans_type,trans_name,i_id]

			#attr_lst = [g_id,t_id,gene_tp,g_name,trans_type,trans_name,int_num,i_id,''.join("level "+str(level))]
			#attr_lst = [g_id,t_id,gene_tp,gene_name,trans_type,trans_name,int_num,i_id,''.join("level "+str(level))]

			#Join attributes to list with ; delimiter
			attr = "; ".join(str(v) for v in attr_lst)

			#Make 2 rows per intron for gene and intron
			VARS1 = [CHR,sorc_b,'gene',START,STOP,SCORE,STR,frame,g_id]
			VARS2 = [CHR,sorc_b,feat_b,START,STOP,SCORE,STR,frame,attr]

			#write to file
			out_GTF.write('\t'.join(VARS1)+';'+'\n'+'\t'.join(VARS2)+';'+'\n')

		elif NAME not in gene_dict:

			continue
			
			##--Make new intron_id--##

			#cont = [str(NAME),str(CHR),str(START),str(STOP),str(STR)]
			#ID = '_'.join(cont)
			##Apply to IDs
			#g_id = "".join("gene_id \""+ID+"\"")
			#g_name="".join("gene_name \""+ID+"\"")
			#t_id = "".join("transcript_id \""+ID+"\"")
			#i_id = "".join("intron_id \""+ID+"\"")
			#trans_name = "".join("transcript_name \""+ID+"\"")

			#attr_lst = [g_id,t_id,'Alyssa_defined',g_name,'Alyssa_defined',trans_name,i_id]

			##Join attributes to list with ; delimiter
			#attr = "; ".join(str(v) for v in attr_lst)

			#Make 2 rows per intron for gene and intron
			#VARS1 = [CHR,sorc_b,'gene',START,STOP,SCORE,STR,'.',g_id]
			#VARS2 = [CHR,sorc_b,feat_b,START,STOP,SCORE,STR,'.',attr]

			#write to file
			#out_GTF.write('\t'.join(VARS1)+';'+'\n'+'\t'.join(VARS2)+';'+'\n')

			#continue

		linenum2+=1
		if linenum2%1000==0: #check status
			status=str(linenum2)+' finished'
			sys.stderr.write(status+'\r')


out_GTF.close()

