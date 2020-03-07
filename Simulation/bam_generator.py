#!/usr/bin/env python
#-*- coding: utf-8 -*-
#===============================================================================
import os
import glob

for i in glob.iglob('*.read1.fastq'):
    name = os.path.basename(i).split('.')[0]
    r1 = i
    r2 = '{}.read2.fastq'.format(name)
    print('bwa mem -M -R \"@RG\\tID:{2}\\tLB:lib1\\tPL:illumina\\tSM:{2}\\tPU:unit1\" -t 8 /path/to/hg19.fa {0} {1} | samtools sort -@ 8 -O BAM -o {2}.bwa.bam - && picard MarkDuplicates I={2}.bwa.bam O={2}.bam M={2}.marked_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true && rm -rf {2}.bwa.bam'.format(r1, r2, name))
