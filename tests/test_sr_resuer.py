"""Tests for scanitd.inference.sr_resuer soft-clipped read rescue workflow."""

import pytest
from pyfaidx import Fasta
from scanitd.base import MappingMode, MicroRegion
from scanitd.inference.sr_resuer import alignment_operation, calculate_variants, update_tdup_ao
from scanitd.ssw import AlignmentMgr


@pytest.fixture()
def sample_fasta(tmp_path):
    """Create a real indexed FASTA file for end-to-end rescue testing."""
    fasta_path = tmp_path / "ref.fa"
    # 0..20: upstream flank, 20..40: duplicated region, 40..60: downstream flank
    seq = "AAAAACCCCCGGGGGTTTTT" + "ACGTACGTACGTACGTACGT" + "TTTTTGGGGGCCCCCTTTTT"
    fasta_path.write_text(f">chr1\n{seq}\n")
    fasta = Fasta(str(fasta_path), sequence_always_upper=True)
    yield fasta
    fasta.close()


class TestCalculateVariants:
    def test_counts_snvs_insertions_and_deletions(self):
        variants = calculate_variants(
            cigar_pair_list=[(3, "M"), (2, "I"), (2, "M"), (1, "D"), (3, "M")],
            reference_seq="ACGTACGGTT",
            query_seq="ACTTTTAGGTT",
            reference_start=0,
            reference_end=9,
            query_start=0,
            query_end=10,
        )
        assert variants["insertions"] == 1
        assert variants["deletions"] == 1
        assert variants["snvs"] == 1  # position 2: G vs T
        assert variants["mismatches"] == 4  # 1 SNV + 2 I + 1 D


class TestAlignmentOperation:
    def test_sm_end_anchored_exact_match_passes(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        ref = "AAAAACCCCCGGGGGTTTTTACGTACGTACGT"
        query = "GGGGGTTTTTACGTACGTACGT"
        assert alignment_operation(mgr, query, ref, MappingMode.SM, mismatches_cutoff=1) is True

    def test_ms_start_anchored_exact_match_passes(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        ref = "ACGTACGTACGTTTTTTGGGGGCCCCCTTTTT"
        query = "ACGTACGTACGTTTTT"
        assert alignment_operation(mgr, query, ref, MappingMode.MS, mismatches_cutoff=1) is True

    def test_unrelated_sequence_is_rejected(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        ref = "AAAAACCCCCGGGGGTTTTTACGTACGTACGT"
        query = "CCCCCCCCCCCCCCCCCCCC"
        assert alignment_operation(mgr, query, ref, MappingMode.SM, mismatches_cutoff=5) is False
        assert alignment_operation(mgr, query, ref, MappingMode.MS, mismatches_cutoff=5) is False

    def test_unanchored_middle_match_is_rejected(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        ref = "AAAAACCCCCGGGGGTTTTTACGTACGTACGT"
        query = "CCCCCGGGGGTTTTT"  # matches middle, not anchored at end or start
        assert alignment_operation(mgr, query, ref, MappingMode.SM, mismatches_cutoff=5) is False
        assert alignment_operation(mgr, query, ref, MappingMode.MS, mismatches_cutoff=5) is False

    def test_excess_mismatches_are_rejected(self):
        mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
        ref = "AAAAACCCCCGGGGGTTTTTACGTACGTACGT"
        # Two substitutions inside the anchored suffix
        query = "GGGGGTTTTTAAGTACGAACGT"
        assert alignment_operation(mgr, query, ref, MappingMode.SM, mismatches_cutoff=1) is False
        assert alignment_operation(mgr, query, ref, MappingMode.SM, mismatches_cutoff=2) is True


class TestUpdateTdupAo:
    def test_blunt_tdup_rescues_sm_and_ms_reads(self, sample_fasta):
        tdup_seq = sample_fasta["chr1"][20:40].seq
        tdup_id = ("chr1", 20, 20, tdup_seq, MicroRegion(""))

        to_be_rescued = {
            ("chr1", 20, MappingMode.SM): [
                "GGGGGTTTTTACGTACGTACGTACGTACGT",  # matches end of ref[0:40]
                "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",  # unrelated, must not rescue
            ],
            ("chr1", 40, MappingMode.MS): [
                "ACGTACGTACGTACGTACGTTTTTTGGGGG",  # matches start of ref[20:60]
            ],
        }

        updated_ao = update_tdup_ao(
            tdup_id=tdup_id,
            original_ao=3,
            genome_fasta=sample_fasta,
            to_be_rescued_sequences=to_be_rescued,
            mismatches_cutoff=1,
        )
        assert updated_ao == 5  # 3 original + 1 SM + 1 MS

    def test_microinsertion_tdup_rescue(self, sample_fasta):
        tdup_seq = sample_fasta["chr1"][20:40].seq
        tdup_id = ("chr1", 20, 20, tdup_seq, MicroRegion("+GATT"))

        to_be_rescued = {
            ("chr1", 20, MappingMode.SM): ["ACGTACGTACGTACGTACGTGATT"],
            ("chr1", 40, MappingMode.MS): ["GATTACGTACGTACGTACGTACGT"],
        }

        updated_ao = update_tdup_ao(
            tdup_id=tdup_id,
            original_ao=2,
            genome_fasta=sample_fasta,
            to_be_rescued_sequences=to_be_rescued,
            mismatches_cutoff=1,
        )
        assert updated_ao == 4

    def test_microhomology_tdup_rescue(self, sample_fasta):
        tdup_seq = sample_fasta["chr1"][20:40].seq
        tdup_id = ("chr1", 20, 20, tdup_seq, MicroRegion("-ACGT"))

        to_be_rescued = {
            # SM uses ref[0:36] (ends with ACGTACGTACGTACGT)
            ("chr1", 20, MappingMode.SM): ["GGGGGTTTTTACGTACGTACGTACGT"],
            # MS uses ref[24:60] (starts with ACGTACGTACGTACGTTTTTT)
            ("chr1", 40, MappingMode.MS): ["ACGTACGTACGTACGTTTTTTGGGGG"],
        }

        updated_ao = update_tdup_ao(
            tdup_id=tdup_id,
            original_ao=1,
            genome_fasta=sample_fasta,
            to_be_rescued_sequences=to_be_rescued,
            mismatches_cutoff=1,
        )
        assert updated_ao == 3
