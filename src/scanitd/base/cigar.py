from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Tuple


@dataclass
class CigarResult:
    cigartuples: List[Tuple[int, int]] = field(default_factory=list)
    cigartuples_without_soft: List[Tuple[int, int]] = field(default_factory=list)
    lt_soft_len: int = 0
    rt_soft_len: int = 0
    read_match: int = 0
    ref_match: int = 0
    indel_len: int = 0
    query_len: int = 0


def parse_cigar(cigar: str) -> CigarResult:
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

    if result.cigartuples[0][0] == 4:  # Left SOFT_CLIP
        result.lt_soft_len = result.cigartuples[0][1]
    if result.cigartuples[-1][0] == 4:  # Right SOFT_CLIP
        result.rt_soft_len = result.cigartuples[-1][1]

    return result
