"""Python ctypes bindings for the vendored Striped Smith-Waterman (SSW) C library.

Adapted from ssw-py (MIT License; see LICENSE-ssw-py.txt).
"""

from __future__ import annotations

import ctypes
from ctypes import POINTER, Structure, byref, c_char_p, c_int, c_int8, c_int32, c_uint8, c_uint16, c_uint32
from pathlib import Path
from typing import ClassVar, NamedTuple

_CIGAR_MAP = "MIDNSHP=X"
_CIGAR_SHIFT = 4
_ALPHABET_SIZE = 5
_SSW_FLAG_ALWAYS_CIGAR = 1

_LIB_PATH = Path(__file__).parent / "ssw.so"
if not _LIB_PATH.is_file():
    _MISSING_LIB_MSG = f"Compiled SSW shared library not found at {_LIB_PATH}; build it with `make -C src/scanitd/ssw`."
    raise ImportError(_MISSING_LIB_MSG)

_lib = ctypes.CDLL(str(_LIB_PATH))


class _SswCtypesResult(Structure):
    _fields_: ClassVar[list[tuple[str, type[ctypes._SimpleCData] | object]]] = [
        ("score1", c_uint16),
        ("score2", c_uint16),
        ("ref_begin1", c_int32),
        ("ref_end1", c_int32),
        ("read_begin1", c_int32),
        ("read_end1", c_int32),
        ("cigar_len", c_int32),
        ("cigar", POINTER(c_uint32)),
    ]


_lib.ssw_align_dna.argtypes = [
    c_char_p,
    c_int32,
    c_char_p,
    c_int32,
    POINTER(c_int8),
    c_int32,
    c_uint8,
    c_uint8,
    c_uint8,
    POINTER(_SswCtypesResult),
]
_lib.ssw_align_dna.restype = c_int

_lib.ssw_free_ctypes_result.argtypes = [POINTER(_SswCtypesResult)]
_lib.ssw_free_ctypes_result.restype = None


class Alignment(NamedTuple):
    """Alignment result returned by :meth:`AlignmentMgr.align`."""

    CIGAR: str
    optimal_score: int
    sub_optimal_score: int
    reference_start: int
    reference_end: int
    read_start: int
    read_end: int
    cigar_pair_list: list[tuple[int, str]]

    def len_cigar(self) -> int:
        """Return the number of CIGAR operations."""
        return len(self.cigar_pair_list)

    def get_cigar_values(self, idx: int) -> tuple[int, str]:
        """Return the ``(length, op)`` CIGAR pair at ``idx``."""
        return self.cigar_pair_list[idx]


def _build_dna_score_matrix(match_score: int, mismatch_penalty: int) -> ctypes.Array[c_int8]:
    """Build a 5x5 int8 DNA substitution matrix for {A, C, G, T, N}."""
    matrix_type = c_int8 * (_ALPHABET_SIZE * _ALPHABET_SIZE)
    values: list[int] = []
    for i in range(_ALPHABET_SIZE):
        for j in range(_ALPHABET_SIZE):
            if i == _ALPHABET_SIZE - 1 or j == _ALPHABET_SIZE - 1:
                values.append(0)
            elif i == j:
                values.append(match_score)
            else:
                values.append(-mismatch_penalty)
    return matrix_type(*values)


class AlignmentMgr:
    """Manager for Striped Smith-Waterman pairwise DNA local alignments."""

    def __init__(
        self,
        match_score: int = 2,
        mismatch_penalty: int = 2,
        gap_open: int = 3,
        gap_extend: int = 1,
    ) -> None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self._score_matrix = _build_dna_score_matrix(match_score, mismatch_penalty)
        self._query: bytes | None = None
        self._reference: bytes | None = None

    def set_read(self, query: str) -> None:
        """Set the query (read) DNA sequence for alignment."""
        if not query:
            msg = "Query sequence must be non-empty."
            raise ValueError(msg)
        self._query = query.encode("ascii")

    def set_reference(self, reference: str) -> None:
        """Set the target (reference) DNA sequence for alignment."""
        if not reference:
            msg = "Reference sequence must be non-empty."
            raise ValueError(msg)
        self._reference = reference.encode("ascii")

    def align(self, gap_open: int | None = None, gap_extension: int | None = None) -> Alignment:
        """Run Striped Smith-Waterman local alignment and return an :class:`Alignment`."""
        if self._query is None or self._reference is None:
            msg = "Read and reference must be set before alignment."
            raise ValueError(msg)

        eff_gap_open = self.gap_open if gap_open is None else gap_open
        eff_gap_ext = self.gap_extend if gap_extension is None else gap_extension
        if eff_gap_open < 0 or eff_gap_ext < 0:
            msg = "Gap penalties must be non-negative integers."
            raise ValueError(msg)

        c_result = _SswCtypesResult()
        try:
            rc = _lib.ssw_align_dna(
                self._query,
                len(self._query),
                self._reference,
                len(self._reference),
                self._score_matrix,
                _ALPHABET_SIZE,
                eff_gap_open,
                eff_gap_ext,
                _SSW_FLAG_ALWAYS_CIGAR,
                byref(c_result),
            )
            if rc != 0:
                msg = f"SSW alignment failed with error code {rc}."
                raise RuntimeError(msg)

            cigar_pair_list: list[tuple[int, str]] = []
            for idx in range(c_result.cigar_len):
                cigar_int = c_result.cigar[idx]
                op_len = int(cigar_int >> _CIGAR_SHIFT)
                op_code = int(cigar_int & 0xF)
                op_char = _CIGAR_MAP[op_code] if op_code < len(_CIGAR_MAP) else "M"
                cigar_pair_list.append((op_len, op_char))

            cigar_str = "".join(f"{op_len}{op_char}" for op_len, op_char in cigar_pair_list)
            return Alignment(
                CIGAR=cigar_str,
                optimal_score=int(c_result.score1),
                sub_optimal_score=int(c_result.score2),
                reference_start=int(c_result.ref_begin1),
                reference_end=int(c_result.ref_end1),
                read_start=int(c_result.read_begin1),
                read_end=int(c_result.read_end1),
                cigar_pair_list=cigar_pair_list,
            )
        finally:
            _lib.ssw_free_ctypes_result(byref(c_result))
