
# ScanITD
[![Build Status](https://travis-ci.org/ylab-hi/ScanITD.svg?branch=master&status=passed)](https://travis-ci.org/ylab-hi/ScanITD)

Introduction
------------
ScanITD: detecting internal tandem duplication with robust variant allele frequency estimation

Getting Started
----------------
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

Prerequisites
----------------
You need Python 3.4 and later to run ScanITD.

### install necessary python packages via anaconda
Install [anaconda](https://www.anaconda.com/download/) (python 3.7) firstly, then install dependent packages via conda in bioconda channel.
```
conda install -c bioconda pysam
conda install -c conda-forge scikit-bio
conda install -c anaconda numpy
conda install -c bioconda samtools  ## samtools/1.0 or newer is required
 ```
Usage
-------------------------
```
ScanITD.py -i input_bam_file -r indexed_refenence_genome_fasta -o output_vcf_filename_prefix [opts]
```

### Options:
```
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        BWA-MEM BAM file
  -r REF, --ref REF     reference genome in FASTA format
  -o OUTPUT, --output OUTPUT
                        output prefix
  -m MAPQ, --mapq MAPQ  minimal MAPQ in BAM for calling ITD (default: 15)
  -c AO, --ao AO        minimal observation count for ITD (default: 4)
  -d DP, --depth DP     minimal depth to call ITD (default: 10)
  -f VAF, --vaf VAF     minimal variant allele frequency (default: 0.1)
  -l ITD_LEN, --len ITD_LEN
                        minimal ITD length to report (default: 10)
  -n MISMATCH           maximum allowed mismatch bases of pairwise local
                        alignment (default: 3)
  -t TARGET, --target TARGET
                        Limit analysis to targets listed in the BED-format
                        file or a samtools region string
  -k, --keep            Kepp the ITD build BAM file
  -v, --version         show program's version number and exit
```  
#### Input:
```	
input_bam_file                    :input WES BAM file. (e.g., wes-seq.bam)
indexed_reference_genome_fasta    :specify reference genome in FASTA format (the reference genome should be indexed)
```
#### Output:
```	
output_vcf_filename_prefix        :specify the prefix of the output vcf file 
```
The name of the output VCF file will be __prefix__.itd.vcf.

The [__Test__](https://github.com/ylab-hi/ScanITD/tree/master/Test "Test directory") directory contains files and scripts for testing the installation.

License
----------------
This project is licensed under [MIT](https://opensource.org/licenses/MIT).

Contact
-----------------
Bug reports or feature requests can be submitted on the [ScanITD Github page](https://github.com/ylab-hi/ScanITD/issues).

Citation
------------------
Wang TY. and Yang R. [ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation](https://doi.org/10.1093/gigascience/giaa089 "ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation").
