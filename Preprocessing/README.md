Prerequisites
----------------
You need Fastqc and Trim Galore to run QC.sh

### install necessary packages via anaconda
Install [anaconda](https://www.anaconda.com/download/) (python 3.7) firstly, then install dependent packages via conda in bioconda channel.
```
conda install -c bioconda trim-galore
conda install -c bioconda fastqc
 ```
 
 Usage
-------------------------
```
trim_galore -o output_folder --paired --fastqc R1.fastq R2.fastq
```
#### Input:
```	
R1.fastq, R2.fastq   :input WES paired-end reads in FASTQ format file.
output_folder        :specify output folder name
```
