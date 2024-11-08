"""Helper functions."""

from __future__ import annotations

import locale
from typing import TYPE_CHECKING, Any

import numpy as np

from scanitd.base import MappingMode, reverse_complement
from scanitd.writer import VCFWriter

if TYPE_CHECKING:
    from pathlib import Path

    from pyfaidx import Fasta
    from pysam import AlignmentFile

    from scanitd.mtype import LoggerType

__all__ = [
    "format_sa_tag",
    "get_insertion_reference_pos",
    "obtain_depth_given_genomic_position",
    "obtain_sa_query_seq_from_ra",
    "parse_target_genomic_coordinates",
    "same_chrom_same_strand_handler",
    "same_chrom_same_strand_mode21_handler",
    "write_events_to_vcf",
]


def same_chrom_same_strand_handler(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    genome_fasta,
    logger,
    microinsertion_cutoff=20,
):
    """Handler for same chrom and same strand."""
    if lt_mode == MappingMode.SM and rt_mode == MappingMode.MS:
        return same_chrom_same_strand_mode21_handler(
            read_lt,
            read_rt,
            lt_mode,
            rt_mode,
            genome_fasta,
            logger,
            microinsertion_cutoff,
        )

    if lt_mode == MappingMode.MS and rt_mode == MappingMode.SM:
        return same_chrom_same_strand_mode21_handler(
            read_rt,
            read_lt,
            rt_mode,
            lt_mode,
            genome_fasta,
            logger,
            microinsertion_cutoff,
            is_reverse=True,
        )
    return None


def same_chrom_same_strand_mode21_handler(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    genome_fasta,
    logger,
    microinsertion_cutoff=20,
    *,
    is_reverse=False,
):
    """Same chrom same strand mode 21 handler."""
    lt_chrm = read_lt.chrom

    if lt_mode == MappingMode.SM and rt_mode == MappingMode.MS:
        target_start = read_rt.ref_start
        target_end = read_lt.ref_end
        target_offset = target_end - target_start

        bp_region_seq_len = read_lt.query_length - read_lt.rt_soft_len - read_rt.lt_soft_len - read_lt.read_match_size - read_rt.read_match_size

        # logger.trace(f"{bp_region_seq_len=}")

        if bp_region_seq_len > microinsertion_cutoff:
            logger.trace(f"{bp_region_seq_len=} > {microinsertion_cutoff=}")
            return None

        query_offset = read_lt.reference_match_size + read_rt.reference_match_size if bp_region_seq_len > 0 else read_lt.reference_match_size + read_rt.reference_match_size + bp_region_seq_len

        lt_bp_seq = obtain_bp_region_seq(
            read_lt,
            lt_mode,
            bp_region_seq_len,
            genome_fasta,
        )
        rt_bp_seq = obtain_bp_region_seq(
            read_rt,
            rt_mode,
            bp_region_seq_len,
            genome_fasta,
        )

        evt_size = query_offset - target_offset

        # logger.trace(f"{evt_size=}, {query_offset=}")
        if evt_size <= 0:  # deletion
            pass
        # reads length < tandem duplication size
        elif evt_size >= query_offset:  # large tandem duplication
            junc_start = read_lt.ref_start
            junc_end = junc_start + evt_size
            if not is_reverse:
                return (
                    "TDUP",
                    (
                        f"{lt_chrm}:{junc_start}",
                        f"{lt_chrm}:{junc_end}",
                        2,
                        1,
                    ),
                    (read_lt.ref_start, read_lt.ref_end),
                    (read_rt.ref_start, read_rt.ref_end),
                    (lt_bp_seq, rt_bp_seq),
                    (read_lt.strand, read_rt.strand),
                )
            # mapping mode and position are not matched here.
            return (
                "TDUP",
                (
                    f"{lt_chrm}:{junc_start}",
                    f"{lt_chrm}:{junc_end}",
                    1,
                    2,
                ),
                (read_rt.ref_start, read_rt.ref_end),
                (read_lt.ref_start, read_lt.ref_end),
                (rt_bp_seq, lt_bp_seq),
                (read_rt.strand, read_lt.strand),
            )
        # read length > tandem duplication size
        else:
            # softclipped length < tandem duplication size (check chimeric read [SM])
            if softclipped_length_and_event_size_checker(
                read_lt,
                lt_mode,
                evt_size,
                bp_region_seq_len,
            ):
                logger.trace("softclipped length < event size: TDUP")
                is_dup = True
            # softclipped length >= tandem duplication size
            # TDUP; Novel Insertion feature: evt_size=0 and bp_region_seq_len>0
            else:
                is_dup = True
                logger.trace("softclipped length >= event size: TDUP")
            if is_dup:
                junc_start = read_lt.ref_start
                junc_end = junc_start + evt_size
                if not is_reverse:
                    return (
                        "TDUP",
                        (
                            f"{lt_chrm}:{junc_start}",
                            f"{lt_chrm}:{junc_end}",
                            2,
                            1,
                        ),
                        (read_lt.ref_start, read_lt.ref_end),
                        (read_rt.ref_start, read_rt.ref_end),
                        (lt_bp_seq, rt_bp_seq),
                        (read_lt.strand, read_rt.strand),
                    )
                # mapping mode and position are not matched here.
                return (
                    "TDUP",
                    (
                        f"{lt_chrm}:{junc_start}",
                        f"{lt_chrm}:{junc_end}",
                        1,
                        2,
                    ),
                    (read_rt.ref_start, read_rt.ref_end),
                    (read_lt.ref_start, read_lt.ref_end),
                    (rt_bp_seq, lt_bp_seq),
                    (read_rt.strand, read_lt.strand),
                )
            return None
    return None


def softclipped_length_and_event_size_checker(
    read,
    mode,
    event_size,
    bp_region_seq_len,
) -> bool:
    """When read length > predicted tandem duplication size.

    check whether the softclipped length is less than the inferred event size

    :param read: a chimeric read
    :param mode: mode for the chimeric read
    :param event_size: event size inferred from 'query_offset - target_offset'

    :type read : Read
    :type mode: int
    :type event_size: int
    :return: whether event_size > softclipped_length (If it is True, it will be a TDUP event)
    :rtype: bool
    """
    return read.lt_soft_len < event_size + bp_region_seq_len if mode == MappingMode.SM else read.rt_soft_len < event_size + bp_region_seq_len


def obtain_bp_region_seq(read, mode, bp_region_seq_len, genome_fasta: Fasta) -> str:
    """Obtain breakpoint region sequence from read.

    :param bp_region_seq_len: the length of the breakpoint region sequence
    :param read:  the chimeirc read
    :param mode: mode for the chimeirc read
    :param genome_fasta: reference genome (pyfaidx.Fasta object)
    :type mode: int
    :return: (putative insertion/microhomology sequence from the read; + means insertion,
        - means microhomology, mode)

    ..note:
        * inserted sequence:
          S-----SM---M    M---MS-----S
          SSSSSXXMMMMM    MMMMMXXSSSSS
        * microhomology:
          S-----SM---M    M---MS-----S
          SSSSSSSXXMMM    MMMXXSSSSSSS
    """
    read_seq = read.query_sequence
    chrom = read.chrom
    bp_region_seq = ""
    # inserted sequence
    if bp_region_seq_len > 0:
        if mode == MappingMode.SM:  # SM
            bp_region_seq = read_seq[: read.lt_soft_len][-bp_region_seq_len:]
        elif mode == MappingMode.MS:  # MS
            bp_region_seq = read_seq[-read.rt_soft_len :][:bp_region_seq_len]
        bp_region_seq = "+" + bp_region_seq
    # microhomology
    elif bp_region_seq_len < 0:
        if mode == MappingMode.SM:  # SM
            bp_region_seq = genome_fasta[chrom][read.ref_start : read.ref_start - bp_region_seq_len].seq

        elif mode == MappingMode.MS:  # MS
            bp_region_seq = genome_fasta[chrom][read.ref_end + bp_region_seq_len : read.ref_end].seq

        bp_region_seq = "-" + bp_region_seq

    return bp_region_seq


def format_sa_tag(in_str: str) -> Any:
    """To keep read.reference_start and start position of SA alignment consistent.

    start position of SA alignment need to subtract 1

    :param in_str: string of supplementary read item in the SA tag
    :return: chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa

    .. note::
         pos_sa, mapq_sa and nm_sa are integral variables now.
    """
    chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = in_str.split(",")
    pos_sa = int(pos_sa) - 1  # type: ignore
    mapq_sa = int(mapq_sa)  # type: ignore
    nm_sa = int(nm_sa)  # type: ignore
    return chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa


def obtain_sa_query_seq_from_ra(
    query_seq_ra: str,
    strand_ra: str,
    strand_sa: str,
) -> str:
    """Helper function to define query_seq for the supplementary alignment.

    :param query_seq_ra: query sequence of representative alignment
    :param strand_ra: direction of representative read (-|+)
    :param strand_sa: direction of supplementary read (-|+)
    :return: query sequence of supplementary alignment
    """
    return query_seq_ra if strand_ra == strand_sa else reverse_complement(query_seq_ra)


def parse_target_genomic_coordinates(input_data):
    """
    Parse genomic coordinates from either BED format or chrom:start-end format.

    Args:
        input_data (str): Either a file path, BED format string, or chrom:start-end format string

    Returns:
        list: List of strings in chrom:start-end format

    Raises:
        ValueError: If input format is invalid
    """

    result = []

    if input_data is None or input_data == "":
        return result

    # Helper function to validate and convert coordinates
    def validate_coordinates(chrom, start, end):
        try:
            start = int(start)
            end = int(end)
            if start < 0 or end < 0 or start >= end:
                raise ValueError
            return f"{chrom}:{start}-{end}"
        except ValueError:
            msg = f"Invalid coordinates: {chrom}:{start}-{end}"
            raise ValueError(msg)

    # Check if input might be a file
    if isinstance(input_data, str):
        # Try to open as file first
        try:
            with open(input_data, encoding=locale.getpreferredencoding(False)) as f:
                lines = f.readlines()
            # Process as BED file
            for line in lines:
                if line.strip() and not line.startswith("#"):
                    fields = line.strip().split("\t")
                    if len(fields) >= 3:
                        result.append(validate_coordinates(fields[0], fields[1], fields[2]))
                    else:
                        msg = f"Invalid BED format in line: {line}"
                        raise ValueError(msg)
            return result
        except FileNotFoundError:
            # Not a file, continue with string processing
            pass

    # Process as string
    entries = input_data.split("\n") if "\n" in input_data else [input_data]

    for entry in entries:
        entry = entry.strip()
        if not entry or entry.startswith("#"):
            continue

        # Check if it's in chrom:start-end format
        if ":" in entry and "-" in entry:
            try:
                chrom, pos = entry.split(":")
                start, end = pos.split("-")
                result.append(validate_coordinates(chrom, start, end))
            except ValueError:
                msg = f"Invalid chrom:start-end format: {entry}"
                raise ValueError(msg)

        # Check if it's in BED format (tab-separated)
        elif "\t" in entry:
            fields = entry.split("\t")
            if len(fields) >= 3:
                result.append(validate_coordinates(fields[0], fields[1], fields[2]))
            else:
                msg = f"Invalid BED format: {entry}"
                raise ValueError(msg)
        else:
            msg = f"Unrecognized format: {entry}"
            raise ValueError(msg)

    return result


def self_loop_checker(insertion_seq: str, left_seq: str, right_seq: str, allowed_mismatched: int = 5) -> tuple[bool, int, str]:
    """Check if the insertion seq is a duplication or not."""

    def count_mismatches_vectorized(str1, str2):
        # Convert strings to numpy arrays of characters
        arr1 = np.array(list(str1))
        arr2 = np.array(list(str2))

        # Get the minimum length to handle strings of different lengths
        min_length = min(len(arr1), len(arr2))

        # Count mismatches in the overlapping part
        mismatches = np.sum(arr1[:min_length] != arr2[:min_length])

        # Add the difference in length as additional mismatches if strings differ in length
        mismatches += abs(len(arr1) - len(arr2))

        return mismatches

    ins_seq = insertion_seq
    ins_len = len(insertion_seq)
    steps = int(ins_len / 2) if ins_len % 2 == 0 else int(ins_len / 2) + 1
    # left rolling
    count = 1
    for _i in range(steps):
        ins_seq = ins_seq[-1:] + ins_seq[0 : len(ins_seq) - 1]
        combo_seq = left_seq[-count:] + right_seq[: ins_len - count]
        count += 1
        if count_mismatches_vectorized(ins_seq, combo_seq) <= allowed_mismatched:
            return True, count, combo_seq
    # right rolling
    ins_seq = insertion_seq
    count = 1
    for _i in range(steps):
        ins_seq = ins_seq[1:] + ins_seq[:1]
        combo_seq = left_seq[-(ins_len - count) :] + right_seq[:count]
        count += 1
        if count_mismatches_vectorized(ins_seq, combo_seq) <= allowed_mismatched:
            return True, len(insertion_seq) - count, combo_seq
    # is a insertion of novel sequence
    return False, 0, ""


def get_insertion_reference_pos(cigar_string, read_pos, insertion_size):
    """
    Find the reference position of an insertion in a CIGAR string.

    Args:
        cigar_string (str): The CIGAR string to parse
        read_pos (int): Starting position of the read on the reference
        insertion_size (int): Size of the insertion to find

    Returns:
        int: Reference position of the insertion, or None if not found
    """
    import re

    # Find the specific insertion pattern in the CIGAR string
    pattern = rf"(\d+)M{insertion_size}I(\d+)M"
    match = re.search(pattern, cigar_string)

    if not match:
        return None

    # Get the number of matches before the insertion
    int(match.group(1))

    # Calculate position in the reference where insertion occurs
    # Start from read_pos and add all operations before the insertion
    current_pos = read_pos

    # Parse the full CIGAR string up to our insertion
    cigar_ops = re.findall(r"(\d+)([MIDNSHPX=])", cigar_string[: match.start() + len(match.group(1)) + 1])

    for length, op in cigar_ops:
        length = int(length)
        # Only these operations consume reference sequence
        if op in "MDN=X":
            current_pos += length

    # Subtract 1 as the insertion occurs between bases
    return current_pos - 1


def obtain_depth_given_genomic_position(
    bam_object: AlignmentFile,
    chrom: str,
    position: int,
    read_mode: MappingMode = MappingMode.SM,
) -> int:
    """Obtain the read depth at genomic position.

    MS: position - 1
    SM: position
    """

    _position = position - 1 if read_mode == MappingMode.MS else position
    return bam_object.count(contig=chrom, start=_position, end=_position + 1)


def write_events_to_vcf(
    output_vcf: Path,
    bam_header: Any,
    events: list,
    logger: LoggerType,
) -> None:
    """Parse splice graph for cliques."""

    vcf_writer = VCFWriter(
        f"{output_vcf}",
        bam_header,
    )
    with vcf_writer.open():
        for idx, event in enumerate(events, 1):
            vcf_writer.write_data(event, f"{idx}")
