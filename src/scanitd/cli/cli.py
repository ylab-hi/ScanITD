"""Console script for scanitd2."""

import typer
import sys
from typing import Optional
from pathlib import Path
from loguru import logger
from enum import Enum
from scanitd import __version__
from scanitd.inference import scan_itd
from scanitd.inference import write_events_to_vcf


def itd_len_type(value: int) -> int:
    """Validate ITD length (assumed from original code)"""
    value = int(value)
    if value < 0:
        raise typer.BadParameter("ITD length must be positive")
    return value


def version_callback(value: bool):
    if value:
        typer.echo(f"ScanITD version: {__version__}")
        raise typer.Exit()


# Define the Enum with the specific log levels
class LogLevel(str, Enum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    DEBUG = "debug"
    TRACE = "trace"


app = typer.Typer(
    context_settings={"help_option_names": ["-h", "--help"]},
    help="ScanITD: detecting internal tandem duplication with robust variant allele frequency estimation",
    add_completion=False,
)


@app.command()
def main(
    input_bam: Path = typer.Option(
        ...,
        "-i",
        "--input",
        help="Aligned BAM file",
        exists=True,
    ),
    ref: Path = typer.Option(
        ...,
        "-r",
        "--ref",
        help="reference genome in FASTA format (with fai index)",
        exists=True,
    ),
    output: str = typer.Option(
        ...,
        "-o",
        "--output",
        help="output VCF file",
    ),
    mapq: int = typer.Option(
        15,
        "-m",
        "--mapq",
        help="minimum MAPQ in BAM for calling ITD",
    ),
    ao: int = typer.Option(
        4,
        "-c",
        "--ao",
        help="minimum observation count for ITD",
    ),
    dp: int = typer.Option(
        10,
        "-d",
        "--depth",
        help="minimum depth to call ITD",
    ),
    vaf: float = typer.Option(
        0.1,
        "-f",
        "--vaf",
        help="minimum variant allele frequency",
    ),
    itd_len: int = typer.Option(
        10,
        "-l",
        "--length",
        callback=itd_len_type,
        help="minimum ITD length to report",
    ),
    mismatch: int = typer.Option(
        3,
        "-n",
        "--mismatches",
        help="maximum allowed mismatches of pairwise local alignment",
    ),
    target: str = typer.Option(
        "",
        "-t",
        "--target",
        help="Limit analysis to targets listed in the BED-format file or a samtools region string",
    ),
    log_level: LogLevel = typer.Option(
        LogLevel.INFO, "-l", "--log-level", help="set the logging level."
    ),
    version: Optional[bool] = typer.Option(
        None,
        "-v",
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Show version and exit",
    ),
):
    """
    Process BAM files to detect internal tandem duplications (ITD).
    Input must be a BWA-MEM aligned BAM file and an indexed reference genome in FASTA format.
    """
    logger.remove()

    logger.add(
        sys.stdout,
        level=log_level.upper(),
        enqueue=True,
        colorize=True,
        backtrace=False,
        diagnose=True,
    )
    event_list, bam_header = scan_itd(
        in_bam_path=input_bam,
        mapq_cutoff=mapq,
        ref_genome=ref,
        target_file=target,
        itd_length_cutoff=itd_len,
        mismatches_cutoff=mismatch,
        logger=logger,
    )

    write_events_to_vcf(output, bam_header, event_list, logger)


if __name__ == "__main__":
    app()
