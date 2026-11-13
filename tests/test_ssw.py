"""Tests for the vendored scanitd.ssw Smith-Waterman alignment package."""

import pytest
from scanitd.ssw import Alignment, AlignmentMgr


class TestAlignmentMgr:
    def test_exact_substring_alignment(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        mgr.set_read("ACGTACGT")
        mgr.set_reference("TTACGTACGTTT")

        res = mgr.align(gap_open=3, gap_extension=1)
        assert isinstance(res, Alignment)
        assert res.CIGAR == "8M"
        assert res.cigar_pair_list == [(8, "M")]
        assert res.optimal_score == 16
        assert res.reference_start == 2
        assert res.reference_end == 9
        assert res.read_start == 0
        assert res.read_end == 7
        assert res.len_cigar() == 1
        assert res.get_cigar_values(0) == (8, "M")

    def test_alignment_with_insertion_in_read(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        mgr.set_read("ACGTTTACGT")
        mgr.set_reference("ACGTACGT")

        res = mgr.align(gap_open=3, gap_extension=1)
        assert res.CIGAR == "3M2I5M"
        assert res.cigar_pair_list == [(3, "M"), (2, "I"), (5, "M")]
        assert res.optimal_score == 12
        assert res.reference_start == 0
        assert res.reference_end == 7
        assert res.read_start == 0
        assert res.read_end == 9

    def test_alignment_with_deletion_in_read(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        mgr.set_read("ACGTACGT")
        mgr.set_reference("ACGTTTACGT")

        res = mgr.align(gap_open=3, gap_extension=1)
        assert res.CIGAR == "3M2D5M"
        assert res.cigar_pair_list == [(3, "M"), (2, "D"), (5, "M")]
        assert res.optimal_score == 12
        assert res.reference_start == 0
        assert res.reference_end == 9
        assert res.read_start == 0
        assert res.read_end == 7

    def test_alignment_with_snv_mismatch_and_n_base(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        mgr.set_read("ACGTNCGT")
        mgr.set_reference("ACGTACGT")

        res = mgr.align(gap_open=3, gap_extension=1)
        assert res.CIGAR == "8M"
        assert res.optimal_score == 14
        assert res.reference_start == 0
        assert res.reference_end == 7
        assert res.read_start == 0
        assert res.read_end == 7

    def test_alignment_with_trailing_softclip_style_mismatch(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        mgr.set_read("ACGTACGTTT")
        mgr.set_reference("ACGTACGT")

        res = mgr.align(gap_open=3, gap_extension=1)
        assert res.CIGAR == "8M"
        assert res.optimal_score == 16
        assert res.reference_start == 0
        assert res.reference_end == 7
        assert res.read_start == 0
        assert res.read_end == 7

    def test_long_sequence_uses_16bit_score_fallback(self):
        seq = "ACGT" * 50  # 200 bp -> max score 400 >= 255
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        mgr.set_read(seq)
        mgr.set_reference("TTTT" + seq + "GGGG")

        res = mgr.align(gap_open=3, gap_extension=1)
        assert res.optimal_score == 400
        assert res.CIGAR == "200M"
        assert res.reference_start == 4
        assert res.reference_end == 203
        assert res.read_start == 0
        assert res.read_end == 199

    def test_no_match_returns_invalid_reference_coords(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        mgr.set_read("AAAA")
        mgr.set_reference("TTTT")

        res = mgr.align(gap_open=3, gap_extension=1)
        assert res.optimal_score == 0
        assert res.reference_start == -1
        assert res.reference_end == -1

    def test_missing_or_empty_inputs_raise_value_error(self):
        mgr = AlignmentMgr()
        with pytest.raises(ValueError, match="Read and reference must be set"):
            mgr.align()

        with pytest.raises(ValueError, match="non-empty"):
            mgr.set_read("")

        with pytest.raises(ValueError, match="non-empty"):
            mgr.set_reference("")

        mgr.set_read("ACGT")
        mgr.set_reference("ACGT")
        with pytest.raises(ValueError, match="non-negative"):
            mgr.align(gap_open=-1)
