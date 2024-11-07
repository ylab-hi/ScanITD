"""Module for writing data to a file."""

from .writer import Writer
from .vcf_writer import VCFWriter

__all__ = ["VCFWriter", "Writer"]
