# Usage Reference

## Command-line interface

```
scanitd [OPTIONS]
```

ScanITD detects internal tandem duplications (ITDs) from a coordinate-sorted BAM
file and writes results in VCF 4.3 format.

---

## Required arguments

| Flag | Short | Type | Description |
|------|-------|------|-------------|
| `--input` | `-i` | PATH | Aligned BAM file (must be indexed) |
| `--ref` | `-r` | PATH | Reference genome FASTA (with .fai index) |
| `--output` | `-o` | TEXT | Output VCF file path |

---

## Optional arguments

### Filtering thresholds

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--mapq` | `-m` | `15` | Minimum MAPQ for a read to be considered |
| `--ao` | `-c` | `4` | Minimum alternate allele observation count |
| `--depth` | `-d` | `10` | Minimum read depth at the locus |
| `--vaf` | `-f` | `0.1` | Minimum variant allele frequency (0–1) |
| `--length` | | `10` | Minimum ITD length in base pairs |

### Alignment sensitivity

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--aln-mismatches` | `-n` | `1` | Max mismatches for soft-read rescue alignment |
| `--ins-mismatches` | | `2` | Max mismatches for insertion-inferred duplication |

### Target region

| Flag | Short | Description |
|------|-------|-------------|
| `--target` | `-t` | Restrict analysis to a BED file or `chr:start-end` region string |

### Other

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--log-level` | `-l` | `info` | Logging verbosity: `trace`, `debug`, `info`, `warning`, `error` |
| `--version` | `-v` | | Print version and exit |
| `--help` | `-h` | | Show help and exit |

---

## Detection strategies

ScanITD uses two complementary strategies to detect ITDs:

### 1. Chimeric-read (SA-tag) detection

Reads carrying an `SA` supplementary alignment tag indicate that a single DNA
fragment mapped to two genomic locations. ScanITD checks whether the two
alignment segments from the same chromosome and strand form a tandem duplication
geometry (one SM + one MS alignment), and if so, records the TDUP coordinates.

### 2. Large-insertion (CIGAR-I) detection

For insertions in the CIGAR string that meet the `--length` threshold, ScanITD
queries the flanking reference sequence and uses a rolling-window comparison
(`self_loop_checker`) to determine whether the inserted sequence is consistent
with a tandem duplication. Events that do not match are classified as novel
insertions (INS).

### 3. Soft-clipped read rescue

Soft-clipped reads that lack an SA tag are aligned against the expected
duplicated reference sequence using Smith-Waterman alignment (ssw-py). Those
whose alignment is end-anchored and passes the mismatch threshold are counted
as additional supporting reads, improving VAF estimates.

---

## Examples

### Basic run

```bash
scanitd -i tumor.bam -r hg38.fa -o tumor_itd.vcf
```

### Sensitive low-VAF detection

```bash
scanitd \
  -i tumor.bam \
  -r hg38.fa \
  -o tumor_itd.vcf \
  --vaf 0.01 \
  --ao 2 \
  --depth 5 \
  --aln-mismatches 2
```

### Targeted region (BED file)

```bash
scanitd -i sample.bam -r hg38.fa -o out.vcf --target panel.bed
```

### Debug mode

```bash
scanitd -i sample.bam -r hg38.fa -o out.vcf --log-level debug
```
