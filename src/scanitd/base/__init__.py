"""Core data structures and utilities for ScanITD.

Exports genomic interval types (:class:`Interval`, :class:`Intervals`),
alignment read representation (:class:`Read`), CIGAR code enums (:class:`CigarCode`),
strand/mode enums (:class:`Strand`, :class:`MappingMode`), and structural variant
event containers (:class:`Event`, :class:`MicroRegion`).
"""

from .basic import (
    CigarCode,
    Event,
    Interval,
    Intervals,
    MappingMode,
    MicroRegion,
    Strand,
)
from .basic_read import Read, reverse_complement

__all__ = [
    "CigarCode",
    "Event",
    "Interval",
    "Intervals",
    "MappingMode",
    "MicroRegion",
    "Read",
    "Strand",
    "reverse_complement",
]
