"""Tests for scanitd.base.basic_read — Read.new() and reverse_complement()."""

from array import array

import pytest

from scanitd.base import MappingMode, Read
from scanitd.base.basic_read import reverse_complement


# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------

class TestReverseComplement:
    def test_simple_uppercase(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        assert reverse_complement("AATT") == "AATT"

    def test_empty_string(self):
        assert reverse_complement("") == ""

    # Bug #13 regressions -----------------------------------------------------

    def test_lowercase_handled(self):
        """Lowercase bases must be complemented correctly, not passed through."""
        assert reverse_complement("atcg") == "cgat"

    def test_mixed_case(self):
        assert reverse_complement("AtCg") == "cGaT"

    def test_n_complemented_to_n(self):
        assert reverse_complement("ANCG") == "CGNT"

    def test_iupac_r_y(self):
        """R (A or G) complements to Y (C or T), and vice versa."""
        assert reverse_complement("R") == "Y"
        assert reverse_complement("Y") == "R"

    def test_iupac_k_m(self):
        assert reverse_complement("K") == "M"
        assert reverse_complement("M") == "K"

    def test_reversal_order(self):
        """Result must be reversed, not just complemented."""
        rc = reverse_complement("AAACCC")
        assert rc == "GGGTTT"


# ---------------------------------------------------------------------------
# Read.new  (classmethod constructor)
# ---------------------------------------------------------------------------

def _qualities(length: int) -> array:
    return array("B", [30] * length)


class TestReadNew:
    def test_simple_match(self):
        seq = "A" * 50
        r = Read.new("read1", "chr1", 1000, "+", "50M", 60, 0, seq, _qualities(50))
        assert r.ref_start == 1000
        assert r.ref_end == 1050       # ref_start + reference_match_size
        assert r.lt_soft_len == 0
        assert r.rt_soft_len == 0
        assert r.query_length == 50

    def test_left_soft_clip(self):
        seq = "A" * 55
        r = Read.new("read2", "chr1", 2000, "+", "5S50M", 60, 0, seq, _qualities(55))
        assert r.lt_soft_len == 5
        assert r.rt_soft_len == 0
        assert r.strand.is_forward()

    def test_right_soft_clip(self):
        seq = "A" * 55
        r = Read.new("read3", "chr1", 3000, "-", "50M5S", 60, 0, seq, _qualities(55))
        assert r.rt_soft_len == 5
        assert r.strand.is_reverse()

    def test_both_soft_clips(self):
        seq = "A" * 60
        r = Read.new("read4", "chr1", 0, "+", "5S50M5S", 60, 0, seq, _qualities(60))
        assert r.lt_soft_len == 5
        assert r.rt_soft_len == 5

    def test_simple_mode_lt_ge_rt_is_sm(self):
        """When lt_soft >= rt_soft the simple_mode is SM."""
        seq = "A" * 60
        r = Read.new("read5", "chr1", 0, "+", "10S50M", 60, 0, seq, _qualities(60))
        assert r.simple_mode == MappingMode.SM

    def test_simple_mode_rt_gt_lt_is_ms(self):
        """When rt_soft > lt_soft the simple_mode is MS."""
        seq = "A" * 60
        r = Read.new("read6", "chr1", 0, "+", "50M10S", 60, 0, seq, _qualities(60))
        assert r.simple_mode == MappingMode.MS

    def test_none_qualities_accepted(self):
        seq = "A" * 50
        r = Read.new("read7", "chr1", 0, "+", "50M", 60, 0, seq, None)
        assert r.query_qualities is None

    def test_hash_and_repr(self):
        seq = "A" * 50
        r = Read.new("read8", "chr1", 0, "+", "50M", 60, 0, seq, None)
        assert isinstance(hash(r), int)
        assert "Read(" in repr(r)

    def test_sms_tuple(self):
        seq = "A" * 60
        r = Read.new("r", "chr1", 0, "+", "5S50M5S", 60, 0, seq, None)
        assert r.sms == (5, 50, 5)

    def test_ref_end_computed_from_cigar(self):
        """ref_end = ref_start + reference_match_size (M+D+N), NOT query length."""
        seq = "A" * 62
        # 50M + 2D + 10M  → ref consumed = 62, read consumed = 60
        r = Read.new("r", "chr1", 100, "+", "50M2D10M", 60, 0, seq, None)
        assert r.ref_end == 100 + 62
