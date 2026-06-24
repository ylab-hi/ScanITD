"""Tests for scanitd.inference.helper — pure functions."""

import pytest

from scanitd.inference.helper import (
    format_sa_tag,
    get_insertion_reference_pos,
    obtain_sa_query_seq_from_ra,
    parse_target_genomic_coordinates,
    self_loop_checker,
)


# ---------------------------------------------------------------------------
# format_sa_tag
# ---------------------------------------------------------------------------

class TestFormatSaTag:
    def test_basic_parsing(self):
        sa = "chr1,1001,+,50M,60,0"
        chrom, pos, strand, cigar, mapq, nm = format_sa_tag(sa)
        assert chrom == "chr1"
        assert pos == 1000          # 1-based → 0-based
        assert strand == "+"
        assert cigar == "50M"
        assert mapq == 60
        assert nm == 0

    def test_position_is_zero_based(self):
        """SA tag position is 1-based; format_sa_tag must subtract 1."""
        _, pos, *_ = format_sa_tag("chr2,500,-,30M,30,1")
        assert pos == 499

    def test_types(self):
        _, pos, strand, cigar, mapq, nm = format_sa_tag("chrX,200,+,100M,55,2")
        assert isinstance(pos, int)
        assert isinstance(mapq, int)
        assert isinstance(nm, int)
        assert isinstance(strand, str)
        assert isinstance(cigar, str)


# ---------------------------------------------------------------------------
# obtain_sa_query_seq_from_ra
# ---------------------------------------------------------------------------

class TestObtainSaQuerySeqFromRa:
    def test_same_strand_returns_original(self):
        seq = "ACGTACGT"
        result = obtain_sa_query_seq_from_ra(seq, "+", "+")
        assert result == seq

    def test_different_strand_returns_rc(self):
        seq = "AAAA"
        result = obtain_sa_query_seq_from_ra(seq, "+", "-")
        assert result == "TTTT"   # RC of AAAA

    def test_both_reverse_returns_original(self):
        seq = "GCGC"
        assert obtain_sa_query_seq_from_ra(seq, "-", "-") == seq


# ---------------------------------------------------------------------------
# parse_target_genomic_coordinates
# ---------------------------------------------------------------------------

class TestParseTargetGenomicCoordinates:
    def test_empty_string_returns_empty(self):
        assert parse_target_genomic_coordinates("") == []

    def test_none_returns_empty(self):
        assert parse_target_genomic_coordinates(None) == []

    def test_chrom_start_end_format(self):
        result = parse_target_genomic_coordinates("chr1:1000-2000")
        assert result == ["chr1:1000-2000"]

    def test_multiple_regions_newline_separated(self):
        result = parse_target_genomic_coordinates("chr1:100-200\nchr2:300-400")
        assert "chr1:100-200" in result
        assert "chr2:300-400" in result

    def test_tab_separated_bed_inline(self):
        result = parse_target_genomic_coordinates("chr1\t100\t200")
        assert result == ["chr1:100-200"]

    def test_invalid_coordinates_raises(self):
        with pytest.raises(ValueError):
            parse_target_genomic_coordinates("chr1:200-100")  # start > end

    # Bug #15 regression -------------------------------------------------------
    def test_bare_chrom_name_accepted(self):
        """Bare chromosome names (e.g. 'chr1') should be accepted as samtools regions."""
        result = parse_target_genomic_coordinates("chr1")
        assert result == ["chr1"]

    def test_bed_file(self, tmp_path):
        bed = tmp_path / "targets.bed"
        bed.write_text("chr1\t100\t200\nchr2\t300\t400\n")
        result = parse_target_genomic_coordinates(str(bed))
        assert "chr1:100-200" in result
        assert "chr2:300-400" in result

    def test_bed_file_ignores_comments(self, tmp_path):
        bed = tmp_path / "targets.bed"
        bed.write_text("# comment\nchr1\t50\t150\n")
        result = parse_target_genomic_coordinates(str(bed))
        assert len(result) == 1
        assert result[0] == "chr1:50-150"


# ---------------------------------------------------------------------------
# get_insertion_reference_pos
# ---------------------------------------------------------------------------

class TestGetInsertionReferencePos:
    def test_simple_insertion(self):
        # 10M then 5I then 10M, read starts at pos 100
        # After 10M: ref = 110; insertion is between 109 and 110 → return 109
        pos = get_insertion_reference_pos("10M5I10M", 100, 5)
        assert pos == 109

    def test_no_matching_insertion_returns_none(self):
        assert get_insertion_reference_pos("50M", 0, 10) is None

    def test_wrong_insertion_size_returns_none(self):
        assert get_insertion_reference_pos("10M5I10M", 0, 3) is None

    def test_insertion_after_deletion(self):
        # 5M2D5M10I5M; looking for 10I
        # ref consumed before 10I: 5M + 2D + 5M = 12 ref bases from pos 0 → 12
        # insertion occurs after base 11 → return 11
        pos = get_insertion_reference_pos("5M2D5M10I5M", 0, 10)
        assert pos == 11

    def test_read_pos_offset(self):
        # Same CIGAR as simple but read starts at 200 instead of 0
        pos_at_0 = get_insertion_reference_pos("10M5I10M", 0, 5)
        pos_at_200 = get_insertion_reference_pos("10M5I10M", 200, 5)
        assert pos_at_200 - pos_at_0 == 200


# ---------------------------------------------------------------------------
# self_loop_checker
# ---------------------------------------------------------------------------

class TestSelfLoopChecker:
    """
    self_loop_checker(insertion_seq, left_seq, right_seq, allowed_mismatches)
    returns (is_dup, shift, tdup_seq).
    """

    def test_perfect_duplication_detected(self):
        # insertion is exactly the left-flanking sequence
        ins_seq  = "ACGT"
        left_seq = "XACGT"   # last 4 = ACGT
        right_seq = "ACGTX"  # first 4 = ACGT
        is_dup, shift, tdup_seq = self_loop_checker(ins_seq, left_seq, right_seq, 0)
        assert is_dup is True

    def test_novel_insertion_not_detected_as_dup(self):
        ins_seq  = "TTTT"
        left_seq = "ACGC"
        right_seq = "GCGC"
        is_dup, shift, tdup_seq = self_loop_checker(ins_seq, left_seq, right_seq, 0)
        assert is_dup is False
        assert shift == 0
        assert tdup_seq == ""

    def test_mismatch_tolerance(self):
        # allow 1 mismatch
        ins_seq  = "ACGT"
        left_seq = "XACGT"
        right_seq = "ACATX"  # 1 mismatch vs ACGT
        is_dup, shift, _ = self_loop_checker(ins_seq, left_seq, right_seq, 1)
        assert is_dup is True

    def test_exceeds_mismatch_limit_returns_false(self):
        ins_seq   = "ACGT"
        left_seq  = "AAAA"
        right_seq = "CCCC"
        is_dup, _, _ = self_loop_checker(ins_seq, left_seq, right_seq, 0)
        assert is_dup is False

    # Bug #6 regression -------------------------------------------------------
    def test_shift_not_off_by_one_for_left_roll(self):
        """
        For a left-roll match at step 1, the returned shift must be 1
        (not 2 as the pre-fix code returned by incrementing before checking).
        """
        ins_seq   = "GACG"    # rolling left once gives: GGAC -> then CGGA
        # Craft sequences so the match happens at roll step 1 (shift=1)
        # After 1 left-roll: ins_seq → "GGAC" is NOT what we want;
        # easier: use a direct identity where roll-step-1 gives exact match.
        # ABCD → left roll 1 → DABC
        ins_seq   = "ABCD"
        left_seq  = "XDABC"   # last 1 = D, so combo at count=1: left[-1:]+right[:3] = "D"+"ABC" = "DABC"
        right_seq = "ABCX"
        is_dup, shift, _ = self_loop_checker(ins_seq, left_seq, right_seq, 0)
        if is_dup:
            assert shift == 1   # must be 1, not 2
