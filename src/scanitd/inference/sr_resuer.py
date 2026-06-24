"""Split-read rescue functions for augmenting tandem duplication allele observations.

Uses Smith-Waterman local alignment (ssw-py) to rescue soft-clipped reads
that support a TDUP event but lack an SA tag, improving allele frequency estimates.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ssw import AlignmentMgr

from scanitd.base import MappingMode

if TYPE_CHECKING:
    from pyfaidx import Fasta

__all__ = [
    "update_tdup_ao",
]


def update_tdup_ao(
    tdup_id: tuple[str, int, int, str],
    original_ao: int,
    genome_fasta: Fasta,
    to_be_rescued_sequences: dict,
    mismatches_cutoff: int = 5,
) -> int:
    """Update the allele observation count for one TDUP event by rescuing soft-clipped reads.

    For both SM and MS breakpoint positions, collects softclipped sequences that
    could not be anchored to a TDUP via an SA tag, aligns them against the
    expected duplicated reference sequence, and counts those passing the mismatch
    threshold as additional supporting reads.

    Args:
        tdup_id: 5-tuple of (chrom, ref_start, tdup_size, tdup_seq, MicroRegion)
            uniquely identifying the TDUP event.
        original_ao: Original SA-tag-derived allele observation count.
        genome_fasta: Reference genome Fasta object for sequence extraction.
        to_be_rescued_sequences: Dict mapping (chrom, position, MappingMode) to
            lists of softclipped sequences to attempt rescue on.
        mismatches_cutoff: Maximum mismatches allowed in a rescue alignment (default: 5).

    Returns:
        int: Updated allele observation count (original_ao + rescued_ao).
    """
    tdup_chrm, tdup_ref_start, tdup_size, _tdup_seq, break_point_region = tdup_id
    tdup_ref_end = tdup_ref_start + tdup_size

    align_mgr = AlignmentMgr(
        match_score=2,
        mismatch_penalty=2,
    )

    rescued_ao = 0
    # SM, tdup_ref_start as breakpoint

    rescued_read_mode = MappingMode.SM
    possible_rescued_id = tdup_chrm, tdup_ref_start, rescued_read_mode

    query_sequences = to_be_rescued_sequences.get(possible_rescued_id)
    if query_sequences:
        __ref_seq_from_genome_partial = genome_fasta[tdup_chrm][tdup_ref_start - tdup_size : tdup_ref_end].seq
        if break_point_region.micro_type == "microinsertion":
            ref_seq_from_genome = __ref_seq_from_genome_partial + break_point_region.sequence
        elif break_point_region.micro_type == "microhomology":
            ref_seq_from_genome = __ref_seq_from_genome_partial[: -break_point_region.length]
        else:
            ref_seq_from_genome = __ref_seq_from_genome_partial

        for query_seq in query_sequences:
            if alignment_operation(
                align_mgr,
                query_seq,
                ref_seq_from_genome,
                rescued_read_mode,
                mismatches_cutoff,
            ):
                rescued_ao += 1

    # MS, tdup_ref_end as breakpoint

    rescued_read_mode = MappingMode.MS
    possible_rescued_id = tdup_chrm, tdup_ref_end, rescued_read_mode

    query_sequences = to_be_rescued_sequences.get(possible_rescued_id)
    if query_sequences:
        __ref_seq_from_genome_partial = genome_fasta[tdup_chrm][tdup_ref_start : tdup_ref_end + tdup_size].seq
        if break_point_region.micro_type == "microinsertion":
            ref_seq_from_genome = break_point_region.sequence + __ref_seq_from_genome_partial
        elif break_point_region.micro_type == "microhomology":
            ref_seq_from_genome = __ref_seq_from_genome_partial[break_point_region.length :]
        else:
            ref_seq_from_genome = __ref_seq_from_genome_partial

        for query_seq in query_sequences:
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
    """Perform Smith-Waterman alignment and check whether the read supports a TDUP breakpoint.

    Aligns query_seq against reference_seq using the provided AlignmentMgr, counts
    mismatches/indels, and verifies that the alignment spans the expected end of the
    sequence (end-anchored for SM mode, start-anchored for MS mode).

    Args:
        align_mgr: Configured ssw.AlignmentMgr instance.
        query_seq: Softclipped read sequence to align.
        reference_seq: Reference sequence spanning the expected duplicated region.
        read_mode: MappingMode indicating which end of the alignment must be anchored.
        mismatches_cutoff: Maximum allowed mismatches (SNVs + indels) (default: 5).

    Returns:
        bool: True if the alignment passes the mismatch filter and is end-anchored correctly.
    """
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

    return (read_mode == MappingMode.SM and reference_end == len(reference_seq) - 1 and query_end == len(query_seq) - 1) or (read_mode == MappingMode.MS and reference_start == 0 and query_start == 0)


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
