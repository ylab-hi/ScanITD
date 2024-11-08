"""VCF Writer class."""

from __future__ import annotations

import datetime
from functools import singledispatchmethod
from typing import IO, TYPE_CHECKING, Any, ClassVar

from loguru import logger

from scanitd import __version__

from .writer import Writer

from scanitd.base import Event


class VCFWriter(Writer):
    """Writer for VCF files.

    .. note::

        1. CHROM: The name of the sequence (typically a chromosome) on which the variation
            is being called. This sequence is usually known as 'the reference sequence',
            i.e. the sequence against which the given sample varies.
        2. POS: The 1-based position of the variation on the given sequence.
        3. ID: The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown
            a ".". Multiple identifiers should be separated by semi-colons without white-space.
        4. REF:The reference base (or bases in the case of an indel) at the given position
            on the given reference sequence.
        5. ALT: The list of alternative alleles at this position.
        6. QUAL: A quality score associated with the inference of the given alleles.
        7. FILTER: A flag indicating which of a given set of filters the variation has
            failed or PASS if all the filters were passed successfully.
        8. INFO: An extensible list of key-value pairs (fields) describing the variation.
            See below for some common fields. Multiple fields are separated by semicolons
            with optional values in the format: <key>=<data>[,data].
        9. FORMAT: An (optional) extensible list of fields for describing the samples.
            See below for some common fields.
        10. SAMPLE: For each (optional) sample described in the file,
            values are given for the fields listed in FORMAT
    """

    num_fields = 10

    reserved_info: ClassVar[dict[str, str]] = {
        "DP": "Integer",
        "OAO": "Integer",
        "AO": "Integer",
        "AF": "Float",
        "SVMETHOD": "String",
        "SVTYPE": "String",
        "SVLEN": "Integer",
        "CHR2": "String",
        "END": "Integer",
        "HOMSEQ": "String",
        "INSSEQ": "String",
        "SEQ": "String",
    }
    reserved_format: ClassVar[dict[str, str]] = {"GT": "String"}
    reserved_alt: ClassVar[list[str]] = [
        "TDUP",
        "INS",
    ]

    description: ClassVar[dict[str, str]] = {
        "DP": "Total read depth at the locus",
        "OAO": "Original alternate allele observations",
        "AO": "Alternate allele observations",
        "AF": "Estimated allele frequency in the range (0,1], representing the ratio of reads showing the alternative allele to all reads",
        "SVTYPE": "The type of event, TDUP, INS.",
        "SVLEN": "Difference in length between REF and ALT alleles",
        "CHR2": "Chromosome for END coordinate in case of a translocation",
        "END": "END coordinate in case of a translocation",
        "SVMETHOD": "Type of approach used to detect SV",
        "GT": "Genotype",
        "INSSEQ": "Sequence of micro-insertion at event breakpoint",
        "HOMSEQ": "Sequence of micro-homology at event breakpoint",
        "SEQ": "Duplication/Insertion sequence",
        "TDUP": "Tandem duplication",
        "INS": "Insertion",
    }

    def __init__(
        self,
        file_path: str,
        bam_header: dict[str, Any],
    ) -> None:
        """Initialize VCFWriter object."""
        super().__init__(file_path)
        self.bam_header = bam_header
        self.sample_name: str = self.file_path.stem
        self.event_id: int = 1

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, fields: list[str], delimiter: str = "\t") -> str:
        """Formatter for writing data."""
        if fields is None or len(fields) != VCFWriter.num_fields:
            logger.warning(
                f"{self.__class__.__name__}: Number of fields is not equal to 10.",
            )
        return delimiter.join(fields) + "\n"

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            logger.warning(f"{self.__class__.__name__}: File is already opened.")
        self.io = self.file_path.open(mode)  # add asyncio support
        if hasattr(self, "write_header"):
            self.write_header()
        return self.io

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
            logger.trace(f"{self.__class__.__name__}: Closing file.")
            self.io.close()  # type: ignore
            self.io = None

    def write_line(self, line: str) -> None:
        """Write line to file."""
        if self.is_opened:
            self.io.write(line)  # type: ignore
        else:
            logger.warning(f"{self.__class__.__name__}: File is not opened.")

    def write_header(self) -> None:
        """Write header to VCF file."""
        self.write_line(self.header)

    @singledispatchmethod
    def write_data(self, data_object: Event, object_id: int) -> None:  # type: ignore
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: Event, event_id: int) -> None:
        """Write Series to VCF file.

        :param data_object: Event to write to file.
        """
        if data_object.event_type is None:
            logger.warning(
                f"{self.__class__.__name__}: No nodes to write to VCF file in Clique {event_id} Series.",
            )
        # _hop_vcf_feature is a dict
        _hop_vcf_feature = get_vcf_features_from_event(
            data_object,
        )
        hop_vcf_feature = vcf_feature_transformer(_hop_vcf_feature, event_id)
        self.write_line(self.formatter(hop_vcf_feature))

    @property
    def header(self) -> str:
        """VCF header provides metadata describing the body of the file."""

        date = datetime.datetime.now().strftime("%Y%m%d")
        source = f"ScanITDv{__version__}"
        reference = f"<CMD={obtain_reference_from_bam_header(self.bam_header)}," 'Description="Alignment parameters">'

        header_lines = [
            "##fileformat=VCFv4.3",
            f"##fileDate={date}",
            f"##source={source}",
            f"##reference={reference}",
        ]
        header_lines += self.get_contigs()

        for _id in VCFWriter.reserved_info:
            _number: str | int = 0 if VCFWriter.reserved_info[_id] == "Flag" else 1
            header_lines.append(
                f"##INFO=<ID={_id},Number={_number},Type={VCFWriter.reserved_info[_id]}," f'Description="{VCFWriter.description[_id]}">',
            )

        for _id in VCFWriter.reserved_format:
            header_lines.append(
                f"##FORMAT=<ID={_id},Number=1,Type={VCFWriter.reserved_format[_id]}," f'Description="{VCFWriter.description[_id]}">',
            )

        for _id in VCFWriter.reserved_alt:
            header_lines.append(
                f'##ALT=<ID={_id},Description="{VCFWriter.description[_id]}">',
            )
        header_lines.append(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}",
        )

        return "\n".join(header_lines) + "\n"

    def get_contigs(self) -> list[str]:
        """Get contigs from BAM file header."""
        return [f"##contig=<ID={contig_dict['SN']},length={contig_dict['LN']}>" for contig_dict in self.bam_header["SQ"]]


def obtain_reference_from_bam_header(bam_header: dict[str, Any]) -> str:
    """Obtain reference info from BAM header."""
    aligners = {
        "CLC",
        "ContextMap2",
        "CRAC",
        "GSNAP",
        "Novoalign",
        "OLego",
        "RUM",
        "Subread",
        "bwa",
        "bowtie",
        "bowtie2",
    }
    avail_aligners = {x.upper() for x in aligners}
    for item in bam_header.get("PG", {}):
        if item["ID"].upper() in avail_aligners:
            return item["CL"]
    return "Unknown"


def get_vcf_features_from_event(
    event: Event,
):
    """Obtain hop vcf features from one series."""

    sv_type = event.event_type
    chrom = event.chrom
    ref_start = event.ref_start
    oao = event.oao
    ao = event.ao
    dp = event.dp
    af = event.af
    event_sequence = event.event_sequence
    end = event.end
    ref_allele = event.ref_allele
    alt_allele = event.alt_allele
    sv_distance = event.event_size
    break_point_region = event.break_point_region

    evt_seq = event_sequence if sv_type in {"TDUP", "INS"} else "."
    micro_insertion = break_point_region.sequence if break_point_region.micro_type == "microinsertion" else "."
    micro_homology = break_point_region.sequence if break_point_region.micro_type == "microhomology" else "."

    return {
        "CHROM": chrom,
        "POS": f"{ref_start + 1}",
        "REF": f"{ref_allele}",
        "ALT": f"{alt_allele}",
        "SVTYPE": sv_type,
        "CHR2": chrom,
        "END": f"{end + 1}",
        "OAO": f"{oao}",
        "AO": f"{ao}",
        "DP": f"{dp}",
        "AF": f"{af:.3g}",
        "SVLEN": f"{sv_distance}",
        "SVMETHOD": "ScanITD2",
        "SEQ": evt_seq,
        "INSSEQ": micro_insertion,
        "HOMSEQ": micro_homology,
    }


def vcf_feature_transformer(feature_dict: dict[str, str], event_id: int) -> list[str]:
    """VCF feature transformer."""
    info_field = (
        f'SVTYPE={feature_dict["SVTYPE"]};OAO={feature_dict["OAO"]};AO={feature_dict["AO"]};'
        f'CHR2={feature_dict["CHR2"]};END={feature_dict["END"]};'
        f'DP={feature_dict["DP"]};AF={feature_dict["AF"]};SVLEN={feature_dict["SVLEN"]};'
        f'INSSEQ={feature_dict["INSSEQ"]};HOMSEQ={feature_dict["HOMSEQ"]};'
        f'SEQ={feature_dict["SEQ"]};SVMETHOD={feature_dict["SVMETHOD"]}'
    )

    return [
        feature_dict["CHROM"],
        feature_dict["POS"],
        str(event_id),
        feature_dict["REF"],
        feature_dict["ALT"],
        ".",
        ".",
        info_field,
        "GT",
        "0/1",
    ]
