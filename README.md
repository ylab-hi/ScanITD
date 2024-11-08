# ScanITD

[![PyPI version](https://img.shields.io/pypi/v/scanitd.svg)](https://pypi.python.org/pypi/scanitd)
[![PyPI - Wheel](https://img.shields.io/pypi/wheel/scanitd)](https://pypi.org/project/scanitd/#files)
[![license](https://img.shields.io/pypi/l/scanitd.svg)](https://github.com/ylab-hi/ScanITD/blob/main/LICENSE)

- **PyPI**

```bash
pip install scanitd
```

- **CONDA** via [Bioconda](https://bioconda.github.io/)

```bash
conda install scanitd
```

# Usage

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

## Credits

- [mit license]: https://opensource.org/licenses/mit
- [pypi]: https://pypi.org/
- [hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python
- [file an issue]: https://github.com/ylab-hi/ScanITD2/issues
- [pip]: https://pip.pypa.io/
- [contributor guide]: CONTRIBUTING.md
- [command-line reference]: https://scanitd.readthedocs.io/en/latest/usage.html
