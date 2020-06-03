Sample data for testing
---
This directory contains one sample dataset for testing purpose to make sure you have installed ScanITD and its dependent packages sucessfully. __test.bam__ is a sliced WXS dataset (in BAM format) contains one FLT3-ITD. 

```
ScanITD.py -i test.bam -r /path/to/hg38.fa -o test 
```
The testing script will generate a VCF file named __test.itd.vcf__
The VCF file contains one ITD (AB=0.12;SVLEN=102)
