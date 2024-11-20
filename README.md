# ScanITD

[![PyPI version](https://img.shields.io/pypi/v/scanitd.svg)](https://pypi.python.org/pypi/scanitd)
[![PyPI - Wheel](https://img.shields.io/pypi/wheel/scanitd)](https://pypi.org/project/scanitd/#files)
[![license](https://img.shields.io/pypi/l/scanitd.svg)](https://github.com/ylab-hi/ScanITD/blob/main/LICENSE)


# 📦 Installation

ScanITD can be installed using pip, the Python package installer.
Follow these steps to install:

1. Ensure you have Python 3.10 or later installed on your system.

2. Create a virtual environment (recommended):

   ```bash
   python -m venv scanitd_env
   source scanitd_env/bin/activate  # On Windows use `scanitd_env\Scripts\activate`
   ```

3. Install ScanITD:

   ```bash
   pip install scanitd
   ```

4. Verify the installation:

   ```bash
   scanitd --help
   ```

# 🛠️ Usage

 Usage: scanitd [OPTIONS]

 ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation
## Required Arguments
* `--input`, `-i` PATH
    - Aligned BAM file
    - Required

* `--ref`, `-r` PATH
    - Reference genome in FASTA format (with fai index)
    - Required

* `--output`, `-o` TEXT
    - Output VCF file
    - Required

## Optional Arguments
* `--mapq`, `-m` INTEGER
    - Minimum MAPQ in BAM for calling ITD
    - Default: 15

* `--ao`, `-c` INTEGER
    - Minimum observation count for ITD
    - Default: 4

* `--depth`, `-d` INTEGER
    - Minimum depth to call ITD
    - Default: 10

* `--vaf`, `-f` FLOAT
    - Minimum variant allele frequency
    - Default: 0.1

* `--length` INTEGER
    - Minimum ITD length to report
    - Default: 10

* `--aln-mismatches`, `-n` INTEGER
    - Maximum allowed mismatches for pairwise local alignment
    - Default: 1

* `--ins-mismatches` INTEGER
    - Maximum allowed mismatches for insertion-inferred duplication
    - Default: 2

* `--target`, `-t` TEXT
    - Limit analysis to targets listed in the BED-format file or a samtools region string

* `--log-level`, `-l` [info|warning|error|debug|trace]
    - Set the logging level
    - Default: info

* `--version`, `-v`
    - Show version and exit

# 📚 Citation

Wang TY. and Yang R. [ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation](https://doi.org/10.1093/gigascience/giaa089 "ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation").
