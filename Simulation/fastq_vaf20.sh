dwgsim -C 80 -1 75 -2 75 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 75bp_100x_WT
dwgsim -C 20 -1 75 -2 75 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 75bp_100x_MT
cat 75bp_100x_WT.bwa.read1.fastq 75bp_100x_MT.bwa.read1.fastq > 75bp_100x.read1.fastq && cat 75bp_100x_WT.bwa.read2.fastq 75bp_100x_MT.bwa.read2.fastq > 75bp_100x.read2.fastq && rm -rf 75bp_100x_WT.bwa.read1.fastq 75bp_100x_MT.bwa.read1.fastq 75bp_100x_WT.bwa.read2.fastq 75bp_100x_MT.bwa.read2.fastq
dwgsim -C 40 -1 75 -2 75 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 75bp_50x_WT
dwgsim -C 10 -1 75 -2 75 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 75bp_50x_MT
cat 75bp_50x_WT.bwa.read1.fastq 75bp_50x_MT.bwa.read1.fastq > 75bp_50x.read1.fastq && cat 75bp_50x_WT.bwa.read2.fastq 75bp_50x_MT.bwa.read2.fastq > 75bp_50x.read2.fastq && rm -rf 75bp_50x_WT.bwa.read1.fastq 75bp_50x_MT.bwa.read1.fastq 75bp_50x_WT.bwa.read2.fastq 75bp_50x_MT.bwa.read2.fastq
dwgsim -C 16 -1 75 -2 75 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 75bp_20x_WT
dwgsim -C 4 -1 75 -2 75 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 75bp_20x_MT
cat 75bp_20x_WT.bwa.read1.fastq 75bp_20x_MT.bwa.read1.fastq > 75bp_20x.read1.fastq && cat 75bp_20x_WT.bwa.read2.fastq 75bp_20x_MT.bwa.read2.fastq > 75bp_20x.read2.fastq && rm -rf 75bp_20x_WT.bwa.read1.fastq 75bp_20x_MT.bwa.read1.fastq 75bp_20x_WT.bwa.read2.fastq 75bp_20x_MT.bwa.read2.fastq
dwgsim -C 80 -1 100 -2 100 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 100bp_100x_WT
dwgsim -C 20 -1 100 -2 100 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 100bp_100x_MT
cat 100bp_100x_WT.bwa.read1.fastq 100bp_100x_MT.bwa.read1.fastq > 100bp_100x.read1.fastq && cat 100bp_100x_WT.bwa.read2.fastq 100bp_100x_MT.bwa.read2.fastq > 100bp_100x.read2.fastq && rm -rf 100bp_100x_WT.bwa.read1.fastq 100bp_100x_MT.bwa.read1.fastq 100bp_100x_WT.bwa.read2.fastq 100bp_100x_MT.bwa.read2.fastq
dwgsim -C 40 -1 100 -2 100 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 100bp_50x_WT
dwgsim -C 10 -1 100 -2 100 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 100bp_50x_MT
cat 100bp_50x_WT.bwa.read1.fastq 100bp_50x_MT.bwa.read1.fastq > 100bp_50x.read1.fastq && cat 100bp_50x_WT.bwa.read2.fastq 100bp_50x_MT.bwa.read2.fastq > 100bp_50x.read2.fastq && rm -rf 100bp_50x_WT.bwa.read1.fastq 100bp_50x_MT.bwa.read1.fastq 100bp_50x_WT.bwa.read2.fastq 100bp_50x_MT.bwa.read2.fastq
dwgsim -C 16 -1 100 -2 100 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 100bp_20x_WT
dwgsim -C 4 -1 100 -2 100 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 100bp_20x_MT
cat 100bp_20x_WT.bwa.read1.fastq 100bp_20x_MT.bwa.read1.fastq > 100bp_20x.read1.fastq && cat 100bp_20x_WT.bwa.read2.fastq 100bp_20x_MT.bwa.read2.fastq > 100bp_20x.read2.fastq && rm -rf 100bp_20x_WT.bwa.read1.fastq 100bp_20x_MT.bwa.read1.fastq 100bp_20x_WT.bwa.read2.fastq 100bp_20x_MT.bwa.read2.fastq
dwgsim -C 80 -1 150 -2 150 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 150bp_100x_WT
dwgsim -C 20 -1 150 -2 150 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 150bp_100x_MT
cat 150bp_100x_WT.bwa.read1.fastq 150bp_100x_MT.bwa.read1.fastq > 150bp_100x.read1.fastq && cat 150bp_100x_WT.bwa.read2.fastq 150bp_100x_MT.bwa.read2.fastq > 150bp_100x.read2.fastq && rm -rf 150bp_100x_WT.bwa.read1.fastq 150bp_100x_MT.bwa.read1.fastq 150bp_100x_WT.bwa.read2.fastq 150bp_100x_MT.bwa.read2.fastq
dwgsim -C 40 -1 150 -2 150 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 150bp_50x_WT
dwgsim -C 10 -1 150 -2 150 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 150bp_50x_MT
cat 150bp_50x_WT.bwa.read1.fastq 150bp_50x_MT.bwa.read1.fastq > 150bp_50x.read1.fastq && cat 150bp_50x_WT.bwa.read2.fastq 150bp_50x_MT.bwa.read2.fastq > 150bp_50x.read2.fastq && rm -rf 150bp_50x_WT.bwa.read1.fastq 150bp_50x_MT.bwa.read1.fastq 150bp_50x_WT.bwa.read2.fastq 150bp_50x_MT.bwa.read2.fastq
dwgsim -C 16 -1 150 -2 150 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 150bp_20x_WT
dwgsim -C 4 -1 150 -2 150 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 150bp_20x_MT
cat 150bp_20x_WT.bwa.read1.fastq 150bp_20x_MT.bwa.read1.fastq > 150bp_20x.read1.fastq && cat 150bp_20x_WT.bwa.read2.fastq 150bp_20x_MT.bwa.read2.fastq > 150bp_20x.read2.fastq && rm -rf 150bp_20x_WT.bwa.read1.fastq 150bp_20x_MT.bwa.read1.fastq 150bp_20x_WT.bwa.read2.fastq 150bp_20x_MT.bwa.read2.fastq
dwgsim -C 80 -1 200 -2 200 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 200bp_100x_WT
dwgsim -C 20 -1 200 -2 200 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 200bp_100x_MT
cat 200bp_100x_WT.bwa.read1.fastq 200bp_100x_MT.bwa.read1.fastq > 200bp_100x.read1.fastq && cat 200bp_100x_WT.bwa.read2.fastq 200bp_100x_MT.bwa.read2.fastq > 200bp_100x.read2.fastq && rm -rf 200bp_100x_WT.bwa.read1.fastq 200bp_100x_MT.bwa.read1.fastq 200bp_100x_WT.bwa.read2.fastq 200bp_100x_MT.bwa.read2.fastq
dwgsim -C 40 -1 200 -2 200 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 200bp_50x_WT
dwgsim -C 10 -1 200 -2 200 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 200bp_50x_MT
cat 200bp_50x_WT.bwa.read1.fastq 200bp_50x_MT.bwa.read1.fastq > 200bp_50x.read1.fastq && cat 200bp_50x_WT.bwa.read2.fastq 200bp_50x_MT.bwa.read2.fastq > 200bp_50x.read2.fastq && rm -rf 200bp_50x_WT.bwa.read1.fastq 200bp_50x_MT.bwa.read1.fastq 200bp_50x_WT.bwa.read2.fastq 200bp_50x_MT.bwa.read2.fastq
dwgsim -C 16 -1 200 -2 200 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.fa 200bp_20x_WT
dwgsim -C 4 -1 200 -2 200 -d 500 -s 50 -c 0 -R 0 -e 0.0005-0.003 -E 0.0005-0.003 chr20.hg19.TDUP.fa 200bp_20x_MT
cat 200bp_20x_WT.bwa.read1.fastq 200bp_20x_MT.bwa.read1.fastq > 200bp_20x.read1.fastq && cat 200bp_20x_WT.bwa.read2.fastq 200bp_20x_MT.bwa.read2.fastq > 200bp_20x.read2.fastq && rm -rf 200bp_20x_WT.bwa.read1.fastq 200bp_20x_MT.bwa.read1.fastq 200bp_20x_WT.bwa.read2.fastq 200bp_20x_MT.bwa.read2.fastq
