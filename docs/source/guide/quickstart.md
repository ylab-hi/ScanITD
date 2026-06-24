# Quick Start

This guide walks through a minimal ScanITD run in under five minutes.

## 1. Prepare your inputs

ScanITD requires three inputs:

| Input | Description |
|-------|-------------|
| Aligned BAM | Coordinate-sorted, indexed BAM from BWA-MEM or compatible aligner |
| Reference FASTA | The genome used for alignment, with `.fai` index |
| Output path | Destination for the VCF results file |

```bash
# Index your BAM (if not already done)
samtools sort sample.bam -o sample.sorted.bam
samtools index sample.sorted.bam

# Index your reference (if not already done)
samtools faidx reference.fa
```

## 2. Run ScanITD

```bash
scanitd \
  --input  sample.sorted.bam \
  --ref    reference.fa \
  --output results.vcf
```

This uses all default thresholds (MAPQ ≥ 15, AO ≥ 4, depth ≥ 10, VAF ≥ 0.1, ITD length ≥ 10 bp).

## 3. Inspect the output

```bash
# Count detected events
grep -v '^#' results.vcf | wc -l

# View first few calls
grep -v '^#' results.vcf | head
```

## Restricting to a target region

To run on a specific gene or genomic region, use `--target`:

```bash
# BED file
scanitd -i sample.bam -r ref.fa -o out.vcf --target targets.bed

# Samtools region string
scanitd -i sample.bam -r ref.fa -o out.vcf --target chr13:28033987-28034454
```

## FLT3-ITD example

FLT3-ITDs are the most clinically relevant ITD class. Focus on the FLT3 locus:

```bash
scanitd \
  --input   sample.bam \
  --ref     hg38.fa \
  --output  flt3_itd.vcf \
  --target  chr13:28033987-28034454 \
  --vaf     0.05 \
  --ao      2
```
