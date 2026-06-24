# Contributing

Thank you for considering contributing to ScanITD!

## Development setup

```bash
git clone https://github.com/ylab-hi/ScanITD.git
cd ScanITD

# Install with dev dependencies using uv
uv sync --dev

# Or with pip
pip install -e ".[dev]"
```

## Running tests

```bash
pytest
```

Run with coverage:

```bash
pytest --cov=scanitd --cov-report=term-missing
```

## Code style

ScanITD uses [ruff](https://docs.astral.sh/ruff/) for linting and formatting.

```bash
# Lint
ruff check src/

# Format
ruff format src/
```

Pre-commit hooks are configured to run these automatically:

```bash
pre-commit install
pre-commit run --all-files
```

## Building documentation

```bash
cd docs
make html
# Open docs/build/html/index.html in your browser
```

## Submitting a pull request

1. Fork the repository and create a branch from `main`.
2. Write tests for any new functionality.
3. Ensure all tests pass and linting is clean.
4. Open a pull request with a clear description of the change.

Please see [CODE\_OF\_CONDUCT.md](https://github.com/ylab-hi/ScanITD/blob/main/CODE_OF_CONDUCT.md) for community standards.
