"""Variant calling pipeline for ScanITD.

Provides the main :func:`scan_itd` entry point for scanning BAM files
and the :func:`write_events_to_vcf` function for VCF output.
"""

from .main import scan_itd
from .helper import write_events_to_vcf

__all__ = ["scan_itd", "write_events_to_vcf"]
