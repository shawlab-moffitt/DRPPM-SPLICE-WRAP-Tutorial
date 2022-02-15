#!/bin/bash


####--------------------------------------------------------####
##    This is a script to generate the reference files and    ##
##          folder structure for the WRAP pipline             ##
####--------------------------------------------------------####


# Example Command: setup.sh [Link to zipped GTF download] [Link to zipped FASTA download] [Prefix] [RemoveCHRflag true/false] [tag: chr or mchr]
# GTF File - Link that will be used to download the GTF with wget
# FASTA link - Link that will be used to download the FASTA with wget
# Prefix - The prefix that will be used when nameing files
# RemoveCHRflag - true/false - Will remove chr from GTF and FASTA and subsequent files
# Tag (OPTIONAL) - Designate what type of chr to remove

# Preferred to start in a folder with a name descriptive of the genome reference files
# Make sure the command line is in an interactive job and the DRPPM software in the PATH

# Download GTF file
#cd GTF
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz
#gunzip gencode.v39lift37.annotation.gtf.gz

#gtf_flag=''
#fasta_flag=''
#prefix=''
#chr_flag='false'
#tag=''

#print_usage() {
#  printf "Usage: ..."
#}

#while getopts 'g:f:p:c:t' flag; do
#  case "${flag}" in
#  	g | --GTFlink) gtf_flag='' ;;
#	f | --FASTAlink) fasta_flag='' ;;
#	p | --prefix) prefix='' ;;
#	c | --CHRflag) chr_flag='false' ;;
#	t | --CHRorMCHR) tag='' ;;
#   h | --help) print_usage
#       exit 1 ;;
#  esac
#done


# Assign gtf file and prefix to be used throughout script
gtf_link=$1
fasta_link=$2
prefix=$3
#prefix=hg19_GRCh37.v39lift37
chrFlag=$4
tag=$5

# Make further sub directories
mkdir GTF
mkdir FASTA
mkdir pipeline
cd GTF
mkdir BED
mkdir CompleteGeneLength
mkdir ExonLength
mkdir ExonNumber

echo "Getting GTF File"

# Download GTF flile
wget ${gtf_link}
gunzip *.gz
gtf_file=$(readlink -f *.gtf)
base=$(basename ${gtf_file})
basenoext=$(basename ${gtf_file} .gtf)

# Remove chr in GTF if flagged
if [ ${chrFlag} = true ]
then
	echo "Removing 'chr' from GTF file"
	drppm -GTFFileAddRemoveChr ${gtf_file} false ${tag} ${basenoext}_nochr.gtf
	gtf_file=$(readlink -f *_nochr.gtf)
	base=$(basename ${gtf_file})
fi

echo "Generating BED files"

# BED files
cd BED
drppm -GTF2BED ${gtf_file} ${prefix}
BED_int=$(readlink -f *.intron.bed)
BED_ex=$(readlink -f *.exon.bed)
BED_gen=$(readlink -f *.gene.bed)

echo "Generating Exon Length Files"

# Exon Length Files
cd ../ExonLength
ln -s ${gtf_file} .
drppm -GTFAnnotateExonLength ${gtf_file} ${prefix}.transcript.exon.length.txt ${prefix}.geneID.exon.length.txt ${prefix}.geneName.exon.length.txt exon.output.txt

echo "Generating Complete Gene Length Files"

# Generate CompleteGeneLength
cd ../CompleteGeneLength
ln -s ${gtf_file} .
drppm -GTFAnnotateGeneLength ${gtf_file} ${prefix}.transcript.gene.length.txt ${prefix}.geneID.gene.length.txt ${prefix}.geneName.gene.length.txt gene.output.txt

echo "Generating Exon Number Files"

# Generate ExonNumber
cd ../ExonNumber
ln -s ${gtf_file} .
drppm -GTFAnnotateNumExon ${gtf_file} ${prefix}.transcript.exon.num.txt ${prefix}.geneID.exon.num.txt ${prefix}.geneName.exon.num.txt ${prefix}.geneName.exonname.num.txt ${prefix}.geneNameExon.num.txt ${prefix}.exon.length.txt

echo "Generating EXON reference Files"

# Generate EXON reference for the htseq pipeline
cd ../../pipeline/
mkdir ExonReference
cd ExonReference
ln -s ${gtf_file} .
drppm -AlternativeJuncGTFFileGenerator ${gtf_file} ${prefix}.exon.gtf ${prefix}.exon.length.txt

echo "Generating Intron reference Files"

# Generate the Intron reference for the htseq pipeline
cd ..
mkdir IntronReference
cd IntronReference
ln -s ${gtf_file} .
ln -s ${BED_int} .
int_dir=$(pwd)
module unload python/2.7.9
module load python/3.9.5
# Make intron GTF
python /share/Lab_Shaw/software/BEDtoGTF.py -b ${BED_int} -g ${gtf_file} -s TimIntron -f intron -o ${int_dir}/${prefix}.intron.gtf
# Make Length Annotation file
python /share/Lab_Shaw/software/GTFtoLengthAnno.py -g ${int_dir}/${prefix}.intron.gtf -o ${int_dir}/${prefix}.intron.length.txt
module unload python/3.9.5
module load python/2.7.9

echo "Generating splicing deficiency Files"

cd ..
mkdir splicing_deficiency
cd splicing_deficiency
echo -e "INTRON_ONLY_BED = ${BED_int}\nEXON_BED = ${BED_ex}\nGENE_BED = ${BED_gen}\nKGXREF = NA\nHG19GTF = ${gtf_file}" > splicing_deficiency.config



# Make STAR Index
module load star/2.6.0c
cd ../../FASTA
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz

echo "Getting FASTA file"

wget ${fasta_link}
gunzip GRCh37.primary_assembly.genome.fa.gz
fasta_file1=$(readlink -f *.genome.fa)
fabasenoext=$(basename ${fasta_file1} .fa)
# Remove chr in GTF if flagged
if [ ${chrFlag} = true ]
then
	echo "Removing 'CHR' from FASTA"
	drppm -FastaAddRemoveChr ${fasta_file1} false ${tag} ${fabasenoext}_nochr.fa
	fasta_file=$(readlink -f *_nochr.fa)
fi

echo "Generating STAR index"

# Perform STAR indexing
mkdir star
cd ..
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir FASTA/star --genomeFastaFiles ${fasta_file} --sjdbGTFfile GTF/${base}