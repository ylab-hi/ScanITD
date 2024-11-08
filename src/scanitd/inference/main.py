#!/usr/bin/env python

import re
from collections import defaultdict
from pathlib import Path

import pysam
from pyfaidx import Fasta, FastaNotFoundError

from scanitd.base import Event, MappingMode, MicroRegion, Read

from .helper import (
    format_sa_tag,
    get_insertion_reference_pos,
    obtain_depth_given_genomic_position,
    obtain_sa_query_seq_from_ra,
    parse_target_genomic_coordinates,
    same_chrom_same_strand_handler,
    self_loop_checker,
)
from .sr_resuer import update_tdup_ao


class BamScanner:
    """BcamScanner scan the bam file and output the result to a file."""

    def __init__(
        self,
        input_bam,
        mapq_cutoff,
        ref_genome,
        microinsertion_cutoff,
        regions,
        logger,
    ) -> None:
        """Initialize the class."""
        self.in_bam_path = input_bam
        self.in_bam_object = pysam.AlignmentFile(input_bam, "rb")

        self.bam_chrom_info = {}
        self.mapq_cutoff = mapq_cutoff
        self.ref_genome = ref_genome.expanduser() if "~" in str(ref_genome) else ref_genome
        self.microinsertion_cutoff = microinsertion_cutoff
        self.regions = regions

        self.logger = logger
        self.header = self._get_bam_header()
        self.total_length = 0

        self.tdup_anchors = {}

        self.genome_fasta = self._get_genome_fasta(self.ref_genome)

    def _check_bam_sort(self, header) -> bool:
        """Check if the bam file is sorted."""
        try:
            return header["HD"]["SO"] == "coordinate"
        except KeyError:
            msg = f"Bam file {self.in_bam_object} is not sorted"
            raise RuntimeError(msg) from KeyError

    def _get_genome_fasta(self, ref_genome) -> Fasta:
        """Get the genome fasta file."""
        try:
            return Fasta(str(ref_genome), sequence_always_upper=True)
        except FastaNotFoundError:
            msg = f"Reference File {ref_genome} is Not Found!"
            raise SystemExit(
                msg,
            ) from FastaNotFoundError

    def _count_chrom_info(self, read):
        """Count the chrom and the chrom start and the chrom end."""
        if read.reference_name in self.bam_chrom_info:
            self.bam_chrom_info[read.reference_name][1] = max(read.reference_end, self.bam_chrom_info[read.reference_name][1])
        else:
            self.bam_chrom_info[read.reference_name] = [
                read.reference_start,
                read.reference_end,
            ]

    def _get_bam_header(self):
        """Get bam header."""
        header = self.in_bam_object.header.as_dict()  # type: ignore
        self._check_bam_sort(header)
        return header

    def iter_bam(self):
        """Iterate the bam file."""
        # supplementary alignment cigarstring extraction
        # key: read.query_name + left S + right S
        self.logger.info("Iter bam file and Extracting primary alignments with SA tags")

        for _region in self.regions:
            for read in self.in_bam_object.fetch(region=_region):
                # XA: Alternative hits https://gist.github.com/crazyhottommy/ed73c7e2daee8383dccb35f224f99714
                if read.has_tag("SA") and read.mapping_quality >= self.mapq_cutoff and not read.is_supplementary and not read.is_secondary and not read.has_tag("XA"):
                    chimeric_aln = read.get_tag("SA")[:-1].split(";")  # type: ignore
                    # skip multi-hop DNA segment
                    if len(chimeric_aln) > 1:
                        continue

                    chrm_ra = read.reference_name
                    pos_ra = read.reference_start
                    strand_ra = "-" if read.is_reverse else "+"

                    cigar_ra = read.cigarstring
                    mapq_ra = read.mapping_quality
                    nm_ra = read.get_tag("NM")
                    seq_ra = read.query_sequence
                    query_qualities_ra = read.query_qualities

                    # only consider the first segment
                    sa_string = chimeric_aln[0]

                    chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = format_sa_tag(sa_string)
                    seq_sa = obtain_sa_query_seq_from_ra(seq_ra, strand_ra, strand_sa)

                    if strand_sa == strand_ra:
                        query_qualities_sa = query_qualities_ra
                    elif query_qualities_ra is None:
                        query_qualities_sa = None
                    else:
                        query_qualities_sa = query_qualities_ra[::-1]

                    if chrm_ra == chrm_sa and strand_ra == strand_sa:
                        read_uno, read_dos = (
                            Read.new(
                                read.query_name,
                                chrm_ra,
                                pos_ra,
                                strand_ra,
                                cigar_ra,
                                mapq_ra,
                                nm_ra,  # type: ignore
                                seq_ra,
                                query_qualities_ra,
                            ),
                            Read.new(
                                read.query_name,
                                chrm_sa,
                                pos_sa,
                                strand_sa,
                                cigar_sa,
                                mapq_sa,
                                nm_sa,  # type: ignore
                                seq_sa,
                                query_qualities_sa,
                            ),
                        )

                        uno_mode = read_uno.simple_mode
                        dos_mode = read_dos.simple_mode

                        event_info = same_chrom_same_strand_handler(
                            read_uno,
                            read_dos,
                            uno_mode,
                            dos_mode,
                            self.genome_fasta,
                            self.logger,
                            self.microinsertion_cutoff,
                        )

                        if event_info is not None:
                            (
                                _event_type,
                                _positions,
                                _read1_info,
                                _read2_info,
                                _insertion_info,
                                _strands,
                            ) = event_info

                            break_point_region = MicroRegion(_insertion_info[0])
                            self.logger.trace(f"{break_point_region=} {read.query_name=}")

                            (
                                tdup_start,
                                tdup_end,
                                _,
                                _,
                            ) = _positions

                            tdup_ref_start = int(tdup_start.split(":")[1])
                            tdup_ref_end = int(tdup_end.split(":")[1])

                            self.tdup_anchors[read.query_name] = (
                                chrm_ra,
                                tdup_ref_start,
                                tdup_ref_end,
                                strand_ra,
                                break_point_region,
                            )

        return self.tdup_anchors


def scan_itd(
    in_bam_path,
    mapq_cutoff,
    ref_genome,
    target_file,
    itd_length_cutoff,
    allowed_mismatches_for_sr_rescue,
    allowed_mismatches_for_insertion,
    logger,
    microinsertion_cutoff: int = 10,
):
    """Main function to run scanbam."""
    regions = parse_target_genomic_coordinates(target_file)

    if len(regions) == 0:
        regions = [None]

    bam_scanner = BamScanner(
        input_bam=Path(in_bam_path),
        mapq_cutoff=mapq_cutoff,
        ref_genome=Path(ref_genome),
        microinsertion_cutoff=microinsertion_cutoff,
        regions=regions,
        logger=logger,
    )
    # iterate over all read of the bam file
    tdup_anchors = bam_scanner.iter_bam()

    bam_object = bam_scanner.in_bam_object
    genome_fasta = bam_scanner.genome_fasta

    to_be_rescued_sequences = defaultdict(list)

    query_reads_total_set = set()

    tdup_ao = defaultdict(int)
    tdup_allele_dict = {}

    ins_ao = defaultdict(int)
    ins_allele_dict = {}

    event_list = []

    for _region in regions:
        try:
            for pileup_column in bam_object.pileup(region=_region, stepper="all", truncate=True):
                # a list of pysam.PileupRead
                for pileup_read in pileup_column.pileups:
                    position_of_pileup_site = pileup_read.query_position
                    read = pileup_read.alignment

                    read_name = read.query_name
                    chrm_ra = read.reference_name
                    pos_ra = read.reference_start
                    strand_ra = "-" if read.is_reverse else "+"

                    cigar_ra = read.cigarstring
                    mapq_ra = read.mapping_quality
                    nm_ra = read.get_tag("NM")
                    seq_ra = read.query_sequence
                    query_qualities_ra = read.query_qualities

                    if position_of_pileup_site and read.mapping_quality >= mapq_cutoff:
                        # in BWA-MEM data, supplmentary alignments will always have H in cigar
                        if "S" in cigar_ra and "H" not in cigar_ra:
                            read_obj = Read.new(
                                read_name,
                                chrm_ra,
                                pos_ra,
                                strand_ra,
                                cigar_ra,
                                mapq_ra,
                                nm_ra,  # type: ignore
                                seq_ra,
                                query_qualities_ra,
                            )
                            # Collect reads with softclipping without TDUP anchors
                            if read_name not in tdup_anchors:
                                read_mode = read_obj.simple_mode
                                if read_mode == MappingMode.MS:
                                    softclipped_sequence = read_obj.query_sequence[-read_obj.rt_soft_len :]
                                    softclipped_position = read_obj.ref_end
                                else:
                                    softclipped_sequence = read_obj.query_sequence[: read_obj.lt_soft_len]
                                    softclipped_position = read_obj.ref_start

                                if read_name not in query_reads_total_set:
                                    to_be_rescued_sequences[(chrm_ra, softclipped_position, read_mode)].append(softclipped_sequence)
                                    query_reads_total_set.add(read_name)

                            # deal with TDUP anchors
                            else:
                                (
                                    _,
                                    tdup_ref_start,
                                    tdup_ref_end,
                                    _,
                                    break_point_region,
                                ) = tdup_anchors[read_name]
                                tdup_size = tdup_ref_end - tdup_ref_start

                                tdup_seq = genome_fasta[chrm_ra][tdup_ref_start:tdup_ref_end].seq

                                tdup_id = (
                                    chrm_ra,
                                    tdup_ref_start,
                                    tdup_size,
                                    tdup_seq,
                                    break_point_region,
                                )

                                ref_allele = genome_fasta[chrm_ra][tdup_ref_start : tdup_ref_start + 1].seq

                                tdup_allele_dict[tdup_id] = (ref_allele, "TDUP")

                                if read_name not in query_reads_total_set:
                                    tdup_ao[tdup_id] += 1
                                    query_reads_total_set.add(read_name)

                        # I in the CIGAR ####
                        if pileup_read.indel >= itd_length_cutoff:
                            insertion_size = pileup_read.indel
                            if re.search(rf"\d+M{insertion_size}I\d+M", cigar_ra):
                                left_seq_from_genome = genome_fasta[chrm_ra][(pileup_column.reference_pos - insertion_size + 2) : pileup_column.reference_pos + 1].seq

                                right_seq_from_genome = genome_fasta[chrm_ra][pileup_column.reference_pos + 1 : (pileup_column.reference_pos + insertion_size)].seq

                                insertion_seq_in_read = seq_ra[position_of_pileup_site + 1 : (position_of_pileup_site + insertion_size + 1)]

                                is_dup, left_shift, tdup_seq = self_loop_checker(
                                    insertion_seq_in_read,
                                    left_seq_from_genome,
                                    right_seq_from_genome,
                                    allowed_mismatches_for_insertion,
                                )

                                if is_dup and get_insertion_reference_pos(cigar_ra, pos_ra, insertion_size) == pileup_column.reference_pos:
                                    tdup_ref_start = pileup_column.reference_pos - left_shift
                                    tdup_id = (
                                        chrm_ra,
                                        tdup_ref_start,
                                        insertion_size,
                                        tdup_seq,
                                        MicroRegion(""),
                                    )

                                    if read_name not in query_reads_total_set:
                                        tdup_ao[tdup_id] += 1
                                        query_reads_total_set.add(read_name)

                                    ref_allele = genome_fasta[chrm_ra][tdup_ref_start : tdup_ref_start + 1].seq
                                    tdup_allele_dict[tdup_id] = (
                                        ref_allele,
                                        "TDUP",
                                    )
                                # Novel sequence insertion
                                else:
                                    ins_id = (
                                        chrm_ra,
                                        pileup_column.reference_pos,
                                        insertion_size,
                                        insertion_seq_in_read,
                                        MicroRegion(""),
                                    )
                                    ref_allele = genome_fasta[chrm_ra][pileup_column.reference_pos : pileup_column.reference_pos + 1].seq
                                    alt_allele = seq_ra[position_of_pileup_site : position_of_pileup_site + insertion_size]
                                    ins_allele_dict[ins_id] = (ref_allele, alt_allele)
                                    ins_ao[ins_id] += 1

        except ValueError as e:
            logger.warning(f"{pileup_column=}, {e=}")
            continue
    for tdup_id in tdup_ao:
        original_ao = tdup_ao[tdup_id]
        new_ao = update_tdup_ao(
            tdup_id,
            original_ao,
            genome_fasta,
            to_be_rescued_sequences,
            allowed_mismatches_for_sr_rescue,
        )

        logger.trace(f"{tdup_id=}, {original_ao=}, {new_ao=}")

        _chrom, _ref_start, _event_size, _event_seq, break_point_region = tdup_id
        depth = obtain_depth_given_genomic_position(bam_object, _chrom, _ref_start)
        ref_allele, alt_allele = tdup_allele_dict[tdup_id]
        event_list.append(Event.new("TDUP", tdup_id, original_ao, new_ao, depth, ref_allele, alt_allele))

    for ins_id in ins_ao:
        _chrom, _ref_start, _event_size, _event_seq, break_point_region = ins_id
        ao = ins_ao[tdup_id]
        depth = obtain_depth_given_genomic_position(bam_object, _chrom, _ref_start)
        ref_allele, alt_allele = ins_allele_dict[ins_id]
        event_list.append(Event.new("INS", ins_id, ao, ao, depth, ref_allele, alt_allele))

    # Sort by chrom, then by reference position
    sorted_event_list = sorted(event_list, key=lambda event: (event.chrom, event.ref_start))
    bam_object.close()

    return sorted_event_list, bam_scanner.header
