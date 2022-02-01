
# DRPPM WRAP Tutorial

## DNA RNA Proteome Post-translational Modifications

[DRPPM](https://github.com/gatechatl/DRPPM) is a library collection of scripts for analyzing DNA/RNA/Proteome/Post-translational Modifications. Below is a example RNA-seq analysis workflow using the DRPPM software on FASTQ files from a study observing the knockdown of the USP7 gene in humans.

## DRPPM Installation

### Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| RSEQC_2.6.4 | STAR_2.7.6a | ucsctools_1.04.00 | GCC_5.5.0 | python_2.7.9 |
| glibc_2.14.1 | samtools_1.1 | fastqc_0.11.7 | R_3.5.1 | Bedtools_2.27.1 |

### Git Checkout

```
git clone git@github.com:gatechatl/DRPPM.git # Clone the repo
```

More information on installation and dependencies can be found on the [DRPPM GitHub Page](https://github.com/gatechatl/DRPPM).

## Pipeline Preparation

### Generate .lst file

1. If starting with FASTQ files, gather a list of FASTQ files and run through -Fastq2FileList function.
   * This step is not needed if starting with BAM files
   * This will pair the FASTQ files on the same line
   * For our example we gathered the path and file name for all the FASTQs in a preliminary [fastq.lst]() file
```
drppm -Fastq2FileList [inputFile] [outputFile]
drppm -Fastq2FileList fastq.lst USP7.lst
```

2. The .lst file will then need to be edited with further information on each line in the format below.
   * Example found here: [USP7.lst]()
   * This was done manually but a script is in development to assist in generation of this file.
```
[Sample_Name]\t[FASTQ_1]\t[FASTQ_2]\t[Read_Length]\t[Forward_or_Reverse] #If using FASTQ files
[Sample_Name]\t[BAM_Files]\t[Read_Length]\t[Forward_or_Reverse]          #If using BAM files
```

### Construct Config File

1. Edit the config file, [hg38_WRAP.config](), which will configure which functions to generate for the pipeline
   * This file contains the location of various files used within the pipeline and some configuration options.
   * Below shows the different arguments in the script. If you want to run the analysis set the boolean to 'false'
```
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
SKIP_HTSEQ_GENE = false           # HTSEQ Gene Quantification
SKIP_JUNCSALVAGER = false         # Splice Junction Calculation
SKIP_RNAEDIT = false              # bam-readcount Function
SKIP_OPTITYPE = false             # HLA Genotyping Prediction
```

### Construct and Run the Step 1 Script

1. An example of this script can be found here: [step1_setup_generate_comprehensive_script.sh]()
2. This script will run the .lst and .config file together and produce further files and folders for the proceeding analysis.
   * The main output script defined in the command consists of a script for each sample containing the pipeline commands.
3. This script contains just one line (seen below) for the -WrappingMyRNAseqAnalysisPipeline fuction
   * You may edit the line in the script and run the script or paste the edited line into the command line to run the function
```
drppm -WrappingMyRNAseqAnalysisPipeline [inputFileLst] [type: FASTQ, BAM] [remapping flag: true or false] [run time config file] [prefix for output folder] [outputShellscript]
drppm -WrappingMyRNAseqAnalysisPipeline USP7.lst FASTQ false hg38_WRAP.config Output execute_everything.sh
```


