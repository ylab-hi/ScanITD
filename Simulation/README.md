Step 1
---
Generate duplication positions using RSVSim
```
Rscript sim.R
```
tandemDuplications.csv will generated in this step.

Step 2
---
Prepare duplication file for svsim
```
prepare_svsim.py
```
tandemDuplications.txt will generated in this step.

Step 3
---
Generate rearranged genome in FASTA format
```
python create_indel_genome.py chr20.hg19.fa tandemDuplications.txt chr20.hg19.TDUP.fa
```
__chr20.hg19.TDUP.fa__ is the rearranged genome generated;
__chr20.hg19.fa__ is the unrranged genome.

Step 4
---
Generate simulated reads in various settings. (reads length, reads depth, variant allele frequency)
```
bash fastq_vaf10.sh # VAF=10%
bash fastq_vaf20.sh # VAF=20%
bash fastq_vaf50.sh # VAF=50%
```

Step 5
---
Map reads using BWA-MEM
```
bam_generator.py
```
Pay attention, you need to replace __/path/to/hg19.fa__ in bam_generator.py using a working hg19.fa
