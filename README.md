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

## Usage

 Usage: scanitd [OPTIONS]

 ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --input           -i      PATH                              Aligned BAM file [default: None] [required]         │
│ *  --ref             -r      PATH                              reference genome in FASTA format (with fai index)   │
│                                                                [default: None]                                     │
│                                                                [required]                                          │
│ *  --output          -o      TEXT                              output VCF file [default: None] [required]          │
│    --mapq            -m      INTEGER                           minimum MAPQ in BAM for calling ITD [default: 15]   │
│    --ao              -c      INTEGER                           minimum observation count for ITD [default: 4]      │
│    --depth           -d      INTEGER                           minimum depth to call ITD [default: 10]             │
│    --vaf             -f      FLOAT                             minimum variant allele frequency [default: 0.1]     │
│    --length                  INTEGER                           minimum ITD length to report [default: 10]          │
│    --aln-mismatches  -n      INTEGER                           maximum allowed mismatches for pairwise local       │
│                                                                alignment                                           │
│                                                                [default: 1]                                        │
│    --ins-mismatches          INTEGER                           maximum allowed mismatches for insertion-inferred   │
│                                                                duplication                                         │
│                                                                [default: 2]                                        │
│    --target          -t      TEXT                              Limit analysis to targets listed in the BED-format  │
│                                                                file or a samtools region string                    │
│    --log-level       -l      [info|warning|error|debug|trace]  set the logging level. [default: info]              │
│    --version         -v                                        Show version and exit                               │
│    --help            -h                                        Show this message and exit.                         │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

## Credits

- [mit license]: https://opensource.org/licenses/mit
- [pypi]: https://pypi.org/
- [hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python
- [file an issue]: https://github.com/ylab-hi/ScanITD2/issues
- [pip]: https://pip.pypa.io/
- [contributor guide]: CONTRIBUTING.md
- [command-line reference]: https://scanitd.readthedocs.io/en/latest/usage.html
