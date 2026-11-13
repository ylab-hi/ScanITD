# Changelog

All notable changes to ScanITD will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.9.3] — 2026-11-13

### Fixed
- Removed the 8,000x `pysam.pileup()` depth cap that caused AF underestimation on deep amplicon data (now uncapped by default via `2**31 - 1`, configurable with `--max-depth`)
- Converted 0-based BED start coordinates (`start + 1`) to 1-based samtools region strings in `parse_target_genomic_coordinates()`, preventing `ValueError: start out of range (-1)` when a BED interval starts at `0`

### Changed
- Vendored the MIT-licensed Striped Smith-Waterman (SSW) C library and `ctypes` Python wrapper under `src/scanitd/ssw`, removing the external `ssw-py` dependency
- Added a Hatchling custom build hook (`src/scanitd/ssw/build_ssw.py`) to compile `ssw.so` from source during wheel builds (with `sse2neon.h` support for ARM64)
- Updated `conda/meta.yaml` for native C extension builds (`{{ compiler('c') }}`, `make`) and removed `ssw-py`

### Added
- Unit tests for `scanitd.ssw.AlignmentMgr` and `Alignment` (`tests/test_ssw.py`)
- Unit and end-to-end split-read rescue tests for `scanitd.inference.sr_resuer` (`tests/test_sr_resuer.py`)

---

## [0.9.2] — 2026-07-16

### Fixed
- Filtering options were not applied (#11)
- Potential error in `inference/main.py` (#10)
- Various bug fixes across the `base`, `inference`, and `writer` modules

### Changed
- Releases now publish to PyPI via GitHub Actions using
  [PyPI trusted publishing](https://docs.pypi.org/trusted-publishers/) (OIDC),
  replacing manual `twine upload`
- Removed `poetry.lock` in favor of a single `pyproject.toml`/`uv.lock`-based
  dependency setup
- Updated Python version requirements and dependencies
- Bumped dev dependencies: `pip` (19.2.3 → 26.1), `wheel` (0.33.6 → 0.38.1)

### Added
- Comprehensive test suite

---

## [0.9.1] — 2024-11-20

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
