# Changelog

All notable changes to ScanITD will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.9.1] — unreleased

### Added
- Full Google-style docstrings across all modules (`base`, `inference`, `writer`, `cli`)
- Sphinx documentation with `autodoc`, `napoleon`, `sphinx-click`, and RTD theme
- User guide pages: Installation, Quick Start, Usage, Output format
- `sphinx-rtd-theme` and `sphinx-click` added to dev dependencies

### Changed
- Dev dependency group reduced from 38 to 19 packages by removing redundant
  flake8 ecosystem (replaced by ruff), `black`, `isort`, `reorder-python-imports`,
  `pynvim`, `cohesion`, `nox`, `xdoctest`, and `safety`
- `[tool.isort]` config section removed from `pyproject.toml`

---

## [0.9.0] — 2024-11-08

### Added
- Initial public release of ScanITD 2.0
- Chimeric-read (SA-tag) based TDUP detection
- Large-insertion (CIGAR-I) based TDUP/INS detection
- Split-read rescue via Smith-Waterman alignment (ssw-py)
- VCF 4.3 output with `OAO`, `AO`, `AF`, `INSSEQ`, `HOMSEQ` INFO fields
- Typer-based CLI with `--target` BED/region support
- Loguru-based structured logging

---

## Citation

Wang TY. and Yang R. [ScanITD: Detecting internal tandem duplication with robust
variant allele frequency estimation](https://doi.org/10.1093/gigascience/giaa089).
*GigaScience*, 2020.
