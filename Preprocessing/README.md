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
ScanITD.py -i input_bam_file -r indexed_refenence_genome_fasta -o output_vcf_filename_prefix [opts]
```
