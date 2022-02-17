
# DRPPM WRAP Tutorial

## DNA RNA Proteome Post-translational Modifications

[DRPPM](https://github.com/gatechatl/DRPPM) is a library collection of scripts for analyzing DNA/RNA/Proteome/Post-translational Modifications. Below is a example RNA-seq analysis workflow using the DRPPM software on FASTQ files from a study observing the knockdown of the USP7 gene in humans.

<p align="center">
  <img src="https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/DRPPM_WRAP_Workflow.PNG" width="600" height="600" />
</p>

## DRPPM Installation

### Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| RSEQC_2.6.4 | STAR_2.7.6a | python_3.9.5 | GCC_5.5.0 | python_2.7.9 |
| samtools_1.1 | fastqc_0.11.7 | Bedtools_2.27.1 | R_3.5.1 | bam-readcount_0.8.0 |

### Git Checkout

```bash
git clone git@github.com:gatechatl/DRPPM.git # Clone the repo
```

More information on installation and dependencies can be found on the [DRPPM GitHub Page](https://github.com/gatechatl/DRPPM).

## Pipeline Preparation

### Generate .lst file

1. If starting with FASTQ files, gather a list of FASTQ files and run through -Fastq2FileList function.
   * This step is not needed if starting with BAM files
   * This will pair the FASTQ files on the same line
   * For our example we gathered the path and file name for all the FASTQs in a preliminary [fastq.lst](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/fastq.lst) file

```bash
drppm -Fastq2FileList [inputFile] [outputFile]
drppm -Fastq2FileList fastq.lst USP7.lst
```

2. The .lst file will then need to be edited with further information on each line in the format below.
   * Example found here: [USP7.lst](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/USP7.lst)
   * This was done manually but a script is in development to assist in generation of this file.

```bash
[Sample_Name]\t[FASTQ_1]\t[FASTQ_2]\t[Read_Length]\t[Forward_or_Reverse] #If using FASTQ files
[Sample_Name]\t[BAM_Files]\t[Read_Length]\t[Forward_or_Reverse]          #If using BAM files
```

### Construct Config File

1. Edit the config file, [hg38_WRAP.config](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/hg38_WRAP.config), which will configure which functions to generate for the pipeline
   * This file contains the location of various files used within the pipeline and some configuration options.
   * Below shows the different arguments in the script. If you want to run the analysis set the boolean to 'false'

```bash
# RSEQC
RSEQC_NOWIG = false               # Adds functions to generate and process BigWig files
## Advance mode ##
SKIP_BAM2FASTQ = false            # Converts BAM Files to FASTQ
SKIP_STAR = false                 # STAR Mapping
SKIP_FASTQC = false               # Generates FASTQC Files
SKIP_RSEQC = false                # RSEQC Functions
SKIP_PSI_PSO_CALC = false         # PSI/PSO Calculation
SKIP_SPLICING_DEFICIENCY = false  # Splicing Deficiency Calculation
SKIP_HTSEQ_EXON_QUANT = false     # HTSEQ Exon Quantification
SKIP_HTSEQ_INTRON_QUANT = false   # HTSEQ Intron Quantification
SKIP_HTSEQ_GENE = false           # HTSEQ Gene Quantification
SKIP_JUNCSALVAGER = false         # Splice Junction Calculation
SKIP_RNAEDIT = false              # bam-readcount Function
SKIP_OPTITYPE = false             # HLA Genotyping Prediction
```

### Advanced Setup for Config File

The current pipeline being demonstrated is referencing the hg38 genome build where all of the reference files within the config file are referencing this build. In anticipation of other genome builds being used, we have generated a script to help setup a reference folder for the genome build of your choice. This script generates GTF and BED files as well as a STAR index directory and other annotation files that are used within the pipeline. Below are steps that should be followed to create this reference directory using a hg19_GRCH37 build example where the 'chr' annotation is not used. The reference files generated below will then need to be input to the main .config file described above.

```bash
# Make a main directory to house your reference files
mkdir hg19_GRCh37.v39lift37_nochr
cd hg19_GRCh37.v39lift37_nochr/
# Copy the setup.sh script into this directory
# Run setup.sh
sh setup.sh [Link to zipped GTF download] [Link to zipped FASTA download] [Prefix] [RemoveCHRflag true/false] [tag: chr or mchr]
sh setup.sh https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz hg19_GRCh37.v39lift37 true chr
``` 

### Construct and Run the Step 1 Script

1. An example of this script can be found here: [step1_setup_generate_comprehensive_script.sh](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/step1_setup_generate_comprehensive_script.sh)
2. This script will run the .lst and .config file together and produce further files and folders for the proceeding analysis.
   * The main output script defined in the command consists of a script for each sample containing the pipeline commands.
3. This script contains just one line (seen below) for the -WrappingMyRNAseqAnalysisPipeline fuction
   * You may edit the line in the script and run the script or paste the edited line into the command line to run the function
4. It is recommended and possibly required to run this script on an high performance computer as it requires a large amount of memory. If you have the option please run in interactive mode.

```bash
drppm -WrappingMyRNAseqAnalysisPipeline [inputFileLst] [type: FASTQ, BAM] [remapping flag: true or false] [run time config file] [prefix for output folder] [outputShellscript]
drppm -WrappingMyRNAseqAnalysisPipeline USP7.lst FASTQ false hg38_WRAP.config Output execute_everything.sh
```

## Running the RNAseq Analysis Pipeline

### Folder Structure

1. The [Step1 Script](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial#construct-and-run-the-step-1-script) generates a variety of .lst files, bash scripts, and folders for each sample that will facilitate the analysis.
2. The 'execute_everything.sh' script that is generated from the [Step1 Script](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial#construct-and-run-the-step-1-script) contains a command to run the pipeline for each sample
3. Below shows what the current folder structure should appear as given one sample:
   * To note: The Intermediate and Output folders contain sub folders for each sample. The pipeline is originally run in the Intermediate folder and the final outputs of the process are transferred to the Output folder for that sample.

```bash
├── alyssa_summary_files.lst
├── execute_everything.sh
├── hg38_WRAP.config
├── htseq_exon_files.lst
├── htseq_gene_files.lst
├── output2matrix.sh
├── psi_pso_files.lst
├── rseqc_files.lst
├── splicing_deficiency_files.lst
├── star_finalout_files.lst
├── step1_setup_generate_comprehensive_script.sh
├── USP7.lst
├── Intermediate
│   ├── JK_A2_1_ERCC_S26
│   │   ├── htseq_exon_level
│   │   │   ├── JK_A2_1_ERCC_S26.exon.htseq.lst
│   │   │   └── JK_A2_1_ERCC_S26.htseq.lst
│   │   ├── htseq_gene_level
│   │   │   └── JK_A2_1_ERCC_S26.htseq.lst
│   │   ├── juncsalvager
│   │   │   └── JK_A2_1_ERCC_S26.juncsalvager.lst
│   │   ├── optitype
│   │   ├── psipso
│   │   │   └── JK_A2_1_ERCC_S26.SJ.file.lst
│   │   ├── qc
│   │   │   └── fastqc
│   │   ├── qc_summary
│   │   ├── rnaediting
│   │   ├── rseqc
│   │   ├── splicingdeficiency
│   │   │   └── JK_A2_1_ERCC_S26_bam_file.lst
│   │   ├── star
│   │       └── JK_A2_1_ERCC_S26_star_file.lst
│   └── JK_A2_1_ERCC_S26.sh
├── Output
│   ├── global_qc_summary
│   │   └── input
│   ├── JK_A2_1_ERCC_S26
│   │   ├── htseq_exon_level
│   │   ├── htseq_gene_level
│   │   ├── juncsalvager
│   │   ├── optitype
│   │   ├── psipso
│   │   ├── qc
│   │   │   └── fastqc
│   │   ├── qc_summary
│   │   ├── rnaediting
│   │   ├── rseqc
│   │   ├── splicingdeficiency
│   │   └── star
```

### Running the Pipeline as a Batch Script

1. Following the setup, the pipeline is run with the 'execute_everything.sh' script.
   * This is regularly done as an array batch job
   * Samples normally take around two days to process but are given a week of wall time to prevent jobs failing due to going over their allotted time.
   * Depending on the system and amount of resources, the lines that set up the parameters of the job may change per user
2. Below is an example input for a batch script used for submission.
   * This example is used to run the pipeline on the Moffitt Cancer Center HPC Colossus Cluster
   * The HPC will be transitioning to a new RED cluster where SLURM will be used to run jobs. Updates will be shown in this README

```bash
#!/bin/bash
#PBS -l walltime=120:00:00                       # 120 hour runtime
#PBS -l nodes=1:ppn=1,pmem=64gb,mem=64gb         # Request 1 node, 1 processor per node, and 64gb of memory
#PBS -t 1-10                                     # Run 10 processes at once

# Path to repositories for software or other required programs
# Possible examples are paths to STAR, DRPPM, RSEQC, and wigToBigWig software or any of the modules listed further down
export PATH=$PATH:PATH/TO/REQUIRED/STAR/DIRECTORY/STAR-2.7.6a/source:PATH/TO/REQUIRED/DRPPM/DIRECTORY/DRPPM-master/export/:PATH/TO/REQUIRED/RSEQC/DIRECTORY/RSeQC-2.6.4/install/share/apps/python-2.7.9/bin/:PATH/TO/REQUIRED/wigToBigWigFunction/DIRECTORY/

# Establish starting directory
# This directory should look like the folder tree shown above
cd /PATH/TO/STARTING/DIRECTORY/

# Assigns an ArrayID (1-10) to each line (p) in the execute_everything.sh script, where each line is a seperate shell script
line=`sed -n "${PBS_ARRAYID}p" execute_everything.sh`

# Essential modules for scripts
# If the user has modularized libraries they may be loaded
# If not, the paths to these programs should be in the PATH section
module load gcc/5.5.0
module load samtools/1.1
module load fastqc/0.11.7
module load python/2.7.9
module load R/3.5.1
module load bedtools2/2.27.1
module load bam-readcount/0.8.0

# Each line/shell script of the excute_everything.sh script if run through here as a seperate but grouped job
${line}
```

## Post-Processing

When the pipeline commences there may be a large number of files output depending on the parameters given on what functions to run. With each sample being run separately from each other, it is importent to gather the data in the end to be able to relate the samples to one another and draw out possible inferences. There are a variety of post-processing script to help in gathering information from all the samples into a single file. The scripts below utilize Python version 3.9.5.

### summarygen.py

The [summarygen.py](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/ExampleRun/software/WRAP_SUMMARY_SCRIPTS/summarygen.py) script takes the input of up to 8 summary output files from a single sample of the pipeline and parses them into a single file. The output file is either a two row or two column file of the variable and the value. When the pipeline commences it will output an output2matrix.sh script that will compbine all the single sample summary files into one matrix. This can be overviewed for quality control of the data to see if there are any outliers or values aspects of the data that should be looked into more carefully.

```python
python summarygen.py -t [TIN] -j [Junction Annotation] -b [BAM Stat] -l [STAR Log] -e [Infer Experiment] -d [Inner Distance] -r [Read Distribution] -n [Intron Summary] -s [Sample Name] -R [Output in Row Format] -C [Output in Column Format]
python summarygen.py -t {SAMPLENAME}.Aligned.sortedByCoord.out.summary.txt -j {SAMPLENAME}_junction_annotation_summary_more.txt -b rseqc_bam_stat_report.txt -l {SAMPLENAME}.STAR.Log.final.out -e {SAMPLENAME}_infer_experiment.txt -d {SAMPLENAME}_inner_distance.txt -r {SAMPLENAME}_read_distribution.txt -n {SAMPLENAME}_intron_summary.txt -s {SAMPLENAME} -R
```

### ReadCountSummary.py

The [ReadCountSummary.py](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/ExampleRun/software/WRAP_SUMMARY_SCRIPTS/ReadCountSummary.py) script takes the input of up to 5 .lst files that contain paths to specific outputs of the pipeline covering HTSEQ gene, exon, and intron counts and FPKM, as well as splicing deficiency and PSI/PSO outputs. The script takes each samples outputs and merges them together based on the feature column to produce a matrix for each type of output file. The HTSEQ and Splicing deficieny outputs will produce files denoting the features as gene symbols and Ensemble IDs. The .lst file contains the location of multiple files per sample, so the output number of files will be greater than the input. 

```python
python ReadCountSummary.py -ht [HTSEQ gene .lst file] -e [HTSEQ exon .lst file] -i [HTSEQ intron .lst file] -s [Splicing Deficiency .lst file] -p [PSI/PSO .lst file]
python ReadCountSummary.py -ht htseq_gene_files.lst -e htseq_exon_files.lst -i htseq_intron_files.lst -s splicing_deficiency_files.lst -p psi_pso_files.lst
```

### ps_PostProcess.py

The [ps_PostProcess.py](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/ExampleRun/software/WRAP_SUMMARY_SCRIPTS/ps_PostProcess.py) script takes the output PSI/PSO files from the ReadCountSummary.py script and filters it for further analysis. This script goes through each line and and observes the values based on a criteria and labels the line as "Keep" or "Remove". The line labeled as "Remove" are either all 'NA', all 1.0, all 0.0, or all 'NA' and 1.0 or all 'NA' and 0.0. This means there is no visible differnce between the samples at the exon making it not as informative. The script outputs three files, the unfiltered (with a column denoting "Keep" or "Remove") and filtered matrix, as well as a filter file which just has two columns of the sample name and the "Keep" of "Remove" column.

```python
python ps_PostProcess.py -p [PSI/PSO Matrix File] -o [Outfile Prefix]
python ps_PostProcess.py -p PSI_Summary.txt -o PSI
```

### htseq_ps_PostProcess.py

The [htseq_ps_PostProcess.py](https://github.com/shawlab-moffitt/DRPPM-WRAP-Tutorial/blob/main/ExampleRun/software/WRAP_SUMMARY_SCRIPTS/htseq_ps_PostProcess.py) script takes the one of the HTSEQ output matrices (using gene symbols) from the ReadCountSummary.py script and one of the PSI/PSO output matrices from either the ReadCountSummary.py or ps_PostProcess.py script. The HTSEQ file will have rows denoted with gene symbols and the PSI/PSO file will denote rows with ExonIDs, so to work with this the gene symbol is extracted from the ExonID in the PSI/PSO file which allows the two files to be merged based on gene symbol. Once matched, the HTSEQ and PSI/PSO samples are separated which produces matrices with the same dimensions that can be compared against each other in further visualizations and analysis. This script outputs three files, an output matrix for bother the HTSEQ and PSI/PSO data as well as a fully merged matrix with the sample columns from the HTSEQ data designated with 'h' and the sample columns from the PSI/PSO data designated with 'p'.

```python
python htseq_ps_PostProcess.py -ht [HTSEQ Summary File] -p [PSI/PSO Summary File] -oh [HTSEQ outfile name] -op [PSI/PSO outfile name] -om [HTSEQ and PSI/PSO merged outfile name]
python htseq_ps_PostProcess.py -ht HTSEQ_Counts_Summary_Symbol.txt -p PSI_Summary_filtered.txt -oh HTSEQcounts_psiProcessed.txt -op PSI_htseqcountsProcessed.txt -om merged_PSI_HTSEQcountsProcessed.txt
```

## Intergation to the DRPPM EASY Shiny App








