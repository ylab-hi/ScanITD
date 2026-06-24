"""Tests for scanitd.base.cigar — parse_cigar()."""

import pytest

from scanitd.base.cigar import CigarResult, parse_cigar


# ---------------------------------------------------------------------------
# Happy-path: well-formed CIGAR strings
# ---------------------------------------------------------------------------

class TestParseCigarHappyPath:
    def test_simple_match(self):
        r = parse_cigar("50M")
        assert r.ref_match == 50
        assert r.read_match == 50
        assert r.query_len == 50
        assert r.lt_soft_len == 0
        assert r.rt_soft_len == 0
        assert r.indel_len == 0

    def test_left_soft_clip(self):
        r = parse_cigar("5S45M")
        assert r.lt_soft_len == 5
        assert r.rt_soft_len == 0
        assert r.query_len == 50   # 5S + 45M

    def test_right_soft_clip(self):
        r = parse_cigar("45M5S")
        assert r.lt_soft_len == 0
        assert r.rt_soft_len == 5

    def test_both_soft_clips(self):
        r = parse_cigar("5S40M5S")
        assert r.lt_soft_len == 5
        assert r.rt_soft_len == 5
        assert r.query_len == 50

    def test_insertion(self):
        r = parse_cigar("10M5I10M")
        # INS subtracts from indel_len
        assert r.indel_len == -5
        assert r.read_match == 25   # 10M + 5I + 10M
        assert r.ref_match == 20    # only M counts for ref

    def test_deletion(self):
        r = parse_cigar("10M5D10M")
        assert r.indel_len == 5
        assert r.ref_match == 25   # 10M + 5D + 10M
        assert r.read_match == 20  # only M counts for read

    def test_complex_cigar(self):
        # from the docstring example: 5S10M2I5M10N10M15S
        r = parse_cigar("5S10M2I5M10N10M15S")
        assert r.lt_soft_len == 5
        assert r.rt_soft_len == 15
        assert r.read_match == 10 + 2 + 5 + 10   # M+I ops
        assert r.ref_match  == 10 + 5 + 10 + 10  # M+N ops (N=ref skip)
        assert r.indel_len  == 10 - 2             # N contrib - I contrib

    def test_cigartuples_without_soft_excludes_softclip(self):
        r = parse_cigar("5S20M5S")
        ops = [op for op, _ in r.cigartuples_without_soft]
        assert 4 not in ops   # op_code 4 == SOFT_CLIP

    def test_returns_cigarresult(self):
        assert isinstance(parse_cigar("10M"), CigarResult)


# ---------------------------------------------------------------------------
# Bug #1 regression: empty / unrecognised CIGAR must not IndexError
# ---------------------------------------------------------------------------

class TestParseCigarEdgeCases:
    def test_empty_string_no_crash(self):
        """parse_cigar('') should return a zeroed CigarResult, not IndexError."""
        r = parse_cigar("")
        assert r.lt_soft_len == 0
        assert r.rt_soft_len == 0
        assert r.cigartuples == []

    def test_star_cigar_no_crash(self):
        """'*' is the SAM unmapped sentinel — no recognised ops, must not crash."""
        r = parse_cigar("*")
        assert r.cigartuples == []
        assert r.lt_soft_len == 0

    def test_hard_clip_not_soft_clip(self):
        """Hard-clip (H) should not set lt_soft_len / rt_soft_len."""
        r = parse_cigar("5H40M5H")
        assert r.lt_soft_len == 0
        assert r.rt_soft_len == 0
