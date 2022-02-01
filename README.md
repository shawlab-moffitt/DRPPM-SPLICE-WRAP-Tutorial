
# DRPPM WRAP Tutorial

## DNA RNA Proteome Post-translational Modifications

[DRPPM](https://github.com/gatechatl/DRPPM) is a library collection of scripts for analyzing DNA/RNA/Proteome/Post-translational Modifications. Below is a example RNA-seq analysis workflow using the DRPPM software on FASTQ files from a study observing the knockdown of the USP7 gene in humans.

## DRPPM Installation

### Git Checkout

```
git clone git@github.com:gatechatl/DRPPM.git # Clone the repo
```

More information on installation and dependencies can be found on the [DRPPM GitHub Page](https://github.com/gatechatl/DRPPM).

## Pipeline Preparation

### Generate .lst file

1. Gather list of FASTQ files and run through -Fastq2FileList function.
   * This will pair the FASTQ files on the same line
```
drppm -Fastq2FileList [inputFile] [outputFile]
```
