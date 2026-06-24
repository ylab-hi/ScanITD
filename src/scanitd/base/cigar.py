"""CIGAR string parsing utilities for BAM alignment records."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Tuple


@dataclass
class CigarResult:
    """Parsed result of a CIGAR string, holding per-operation statistics."""

    cigartuples: List[Tuple[int, int]] = field(default_factory=list)
    #: All (op_code, length) CIGAR operation pairs.
    cigartuples_without_soft: List[Tuple[int, int]] = field(default_factory=list)
    #: CIGAR pairs excluding soft-clip operations.
    lt_soft_len: int = 0
    #: Length of left soft-clipped bases.
    rt_soft_len: int = 0
    #: Length of right soft-clipped bases.
    read_match: int = 0
    #: Bases consumed from the read (M+I).
    ref_match: int = 0
    #: Bases consumed from the reference (M+D+N).
    indel_len: int = 0
    #: Net indel length relative to reference (D+N-I).
    query_len: int = 0
    #: Total query length including soft-clips (M+I+S).


def parse_cigar(cigar: str) -> CigarResult:
    """Parse a SAM CIGAR string into a CigarResult.

    Args:
        cigar: CIGAR string (e.g. ``5S10M2I5M``).

    Returns:
        A CigarResult populated with per-operation counts and soft-clip lengths.
    """
    result = CigarResult()
    cigar_pattern = re.compile(r"(\d+)([MIDNSHP=XB])")
    cigar_tuples = cigar_pattern.findall(cigar)

    for op_len, op in cigar_tuples:
        op_len = int(op_len)
        op_code = "MIDNSHP=XB".index(op)
        result.cigartuples.append((op_code, op_len))

        if op_code == 0:  # MATCH
            result.ref_match += op_len
            result.read_match += op_len
            result.query_len += op_len
            result.cigartuples_without_soft.append((op_code, op_len))
        elif op_code == 1:  # INS
            result.indel_len -= op_len
            result.read_match += op_len
            result.query_len += op_len
            result.cigartuples_without_soft.append((op_code, op_len))
        elif op_code in (2, 3):  # DEL or REF_SKIP
            result.indel_len += op_len
            result.ref_match += op_len
            result.cigartuples_without_soft.append((op_code, op_len))
        elif op_code == 4:  # SOFT_CLIP
            result.query_len += op_len

    if result.cigartuples and result.cigartuples[0][0] == 4:  # Left SOFT_CLIP
        result.lt_soft_len = result.cigartuples[0][1]
    if result.cigartuples and result.cigartuples[-1][0] == 4:  # Right SOFT_CLIP
        result.rt_soft_len = result.cigartuples[-1][1]

    return result
