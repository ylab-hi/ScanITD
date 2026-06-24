# Installation

ScanITD requires Python 3.9 or later and works on Linux and macOS.

## Using pip

```bash
pip install scanitd
```

## Using uv (recommended)

```bash
uv tool install scanitd
```

## From source

```bash
git clone https://github.com/ylab-hi/ScanITD.git
cd ScanITD
pip install -e .
```

## Prerequisites

| Requirement | Minimum version | Notes |
|-------------|----------------|-------|
| Python      | 3.9             |       |
| pysam       | 0.22.0          | Requires htslib |
| samtools    | —               | BAM must be coordinate-sorted & indexed |

Your BAM file must be **coordinate-sorted** and accompanied by a `.bai` index.
Your reference FASTA must have a `.fai` index (generate with `samtools faidx ref.fa`).

## Verifying the installation

```bash
scanitd --help
```

You should see the help text listing all available options.
