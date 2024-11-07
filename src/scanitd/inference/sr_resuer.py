"""SR rescuer functions."""

from __future__ import annotations
from ssw import AlignmentMgr

from scanitd.base import MappingMode


__all__ = [
    "update_tdup_ao",
]


def update_tdup_ao(
    tdup_id: tuple[str, int, int, str],
    original_ao: int,
    tdup_anchors_sequences: dict,
    to_be_rescued_sequences: dict,
    mismatches_cutoff: int = 5,
) -> int:
    """Update AO for one TDUP event."""
    tdup_chrm, tdup_ref_start, tdup_size, tdup_seq = tdup_id
    tdup_ref_end = tdup_ref_start + tdup_size

    left_softclipped_sequence, right_softclipped_sequence = tdup_anchors_sequences[
        tdup_id
    ]

    align_mgr = AlignmentMgr(
        match_score=2,
        mismatch_penalty=2,
    )

    ref_seq_from_genome = tdup_seq

    rescued_ao = 0
    for rescued_id in to_be_rescued_sequences:
        chrm_rescued, rescued_softclipped_position, rescued_read_mode = rescued_id
        if chrm_rescued == tdup_chrm:
            ref_seq_from_read = ""
            if (
                rescued_read_mode == MappingMode.SM
                and rescued_softclipped_position == tdup_ref_start
            ):
                ref_seq_from_read = left_softclipped_sequence
            elif (
                rescued_read_mode == MappingMode.MS
                and rescued_softclipped_position == tdup_ref_end
            ):
                ref_seq_from_read = right_softclipped_sequence

            if ref_seq_from_read:
                for query_seq in to_be_rescued_sequences[rescued_id]:
                    if alignment_operation(
                        align_mgr,
                        query_seq,
                        ref_seq_from_genome,
                        rescued_read_mode,
                        mismatches_cutoff,
                    ):
                        rescued_ao += 1

    return original_ao + rescued_ao


def alignment_operation(
    align_mgr: AlignmentMgr,
    query_seq: str,
    reference_seq: str,
    read_mode: MappingMode,
    mismatches_cutoff: int = 5,
):
    align_mgr.set_read(query_seq)
    align_mgr.set_reference(reference_seq)

    # Compute the alignment
    alignment = align_mgr.align(
        gap_open=3,
        gap_extension=1,
    )

    reference_start = alignment.reference_start
    reference_end = alignment.reference_end
    query_start = alignment.read_start
    query_end = alignment.read_end
    variants = calculate_variants(
        alignment.cigar_pair_list,
        reference_seq,
        query_seq,
        reference_start,
        reference_end,
        query_start,
        query_end,
    )

    if variants["mismatches"] > mismatches_cutoff:
        return False

    return (
        read_mode == MappingMode.SM
        and reference_end == len(reference_seq) - 1
        and query_end == len(query_seq) - 1
    ) or (read_mode == MappingMode.MS and reference_start == 0 and query_start == 0)


def calculate_variants(
    cigar_pair_list: list,
    reference_seq: str,
    query_seq: str,
    reference_start: int,
    reference_end: int,
    query_start: int,
    query_end: int,
):
    """
    Calculate the number of indels, deletions, and SNVs from CIGAR string and sequences.

    Parameters:
        cigar_string (str): CIGAR string (e.g., "9M1I5M")
        read_seq (str): The read sequence
        ref_seq (str): The reference sequence
        ref_start (int): Start position in reference (0-based)
        ref_end (int): End position in reference (0-based)
        read_start (int): Start position in read (0-based)
        read_end (int): End position in read (0-based)

    Returns:
        dict: Dictionary containing counts of indels, deletions, and SNVs, and details
    """
    # Initialize counters and details
    variants = {
        "insertions": 0,  # Only insertions
        "deletions": 0,  # Only deletions
        "snvs": 0,  # Single nucleotide variants
        "mismatches": 0,  # mismatches bases including snvs, deletions, and insertions
    }

    # Initialize positions
    read_pos = query_start
    ref_pos = reference_start

    for length, op in cigar_pair_list:
        if op == "M":  # Match/mismatch
            # Compare sequences to find SNVs
            for i in range(length):
                if read_pos + i < query_end and ref_pos + i < reference_end:
                    read_base = query_seq[read_pos + i]
                    ref_base = reference_seq[ref_pos + i]
                    if read_base != ref_base:
                        variants["snvs"] += 1
                        variants["mismatches"] += 1
            read_pos += length
            ref_pos += length

        elif op == "I":  # Insertion
            variants["insertions"] += 1
            variants["mismatches"] += length
            read_pos += length

        elif op == "D":  # Deletion
            variants["deletions"] += 1
            variants["mismatches"] += length
            ref_pos += length

    return variants
