"""Integration regression tests for scanitd.inference.main."""

import pysam
import pytest
from loguru import logger

from scanitd.inference.main import scan_itd


def _build_deep_itd_bam(tmp_path, total_reads: int = 8500, alt_reads: int = 850):
    """Create a synthetic indexed FASTA, BAM (>8,000x depth), and start=0 BED file."""
    # Reference: 200 bp non-repetitive sequence around position 100
    ref_seq = (
        "ACGTACGTACGTACGTACGT" * 4  # 0..80
        + "GCTTAGCTAGGCTTAGCTAG"    # 80..100
        + "TCGATCGATCGATCGATCGA"    # 100..120
        + "TTAACCGGTTAACCGGTTAA" * 4  # 120..200
    )
    fasta_path = tmp_path / "ref.fa"
    fasta_path.write_text(f">chr1\n{ref_seq}\n")
    pysam.faidx(str(fasta_path))

    # Duplicated 15-bp unit right before pos 100: ref_seq[85:100]
    itd_unit = ref_seq[85:100]
    ref_read_seq = ref_seq[60:140]  # 80M spanning 60..140
    alt_read_seq = ref_seq[60:100] + itd_unit + ref_seq[100:140]  # 40M15I40M

    bam_path = tmp_path / "deep.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": len(ref_seq)}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
        for i in range(total_reads):
            is_alt = i < alt_reads
            seg = pysam.AlignedSegment()
            seg.query_name = f"read_{i}"
            seg.query_sequence = alt_read_seq if is_alt else ref_read_seq
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = 60
            seg.mapping_quality = 60
            seg.cigarstring = "40M15I40M" if is_alt else "80M"
            seg.set_tag("NM", len(itd_unit) if is_alt else 0)
            out_bam.write(seg)

    pysam.index(str(bam_path))

    bed_path = tmp_path / "target.bed"
    bed_path.write_text("chr1\t0\t200\n")

    return bam_path, fasta_path, bed_path


class TestScanItdDepthAndBed:
    def test_deep_pileup_not_capped_at_8000_and_bed_start_zero_works(self, tmp_path):
        bam_path, fasta_path, bed_path = _build_deep_itd_bam(tmp_path, total_reads=8500, alt_reads=850)

        events, _ = scan_itd(
            in_bam_path=bam_path,
            mapq_cutoff=15,
            ref_genome=fasta_path,
            target_file=str(bed_path),
            itd_length_cutoff=10,
            allowed_mismatches_for_sr_rescue=1,
            allowed_mismatches_for_insertion=2,
            logger=logger,
            max_depth=0,
        )

        assert len(events) == 1
        ev = events[0]
        assert ev.event_type == "TDUP"
        assert ev.dp == 8500
        assert ev.ao == 850
        assert ev.af == pytest.approx(0.10, abs=1e-6)

    def test_explicit_max_depth_caps_pileup_reads(self, tmp_path):
        bam_path, fasta_path, bed_path = _build_deep_itd_bam(tmp_path, total_reads=8500, alt_reads=850)

        # Write alt reads at the end by passing alt_reads=500 at start and capping at 400
        # All 850 alt reads are written first (i < 850), so capping at 400 caps AO at 400.
        events, _ = scan_itd(
            in_bam_path=bam_path,
            mapq_cutoff=15,
            ref_genome=fasta_path,
            target_file=str(bed_path),
            itd_length_cutoff=10,
            allowed_mismatches_for_sr_rescue=1,
            allowed_mismatches_for_insertion=2,
            logger=logger,
            max_depth=400,
        )

        assert len(events) == 1
        assert events[0].ao == 400
        assert events[0].dp == 8500

    def test_negative_max_depth_raises(self, tmp_path):
        bam_path, fasta_path, bed_path = _build_deep_itd_bam(tmp_path, total_reads=20, alt_reads=5)
        with pytest.raises(ValueError, match="max_depth must be non-negative"):
            scan_itd(
                in_bam_path=bam_path,
                mapq_cutoff=15,
                ref_genome=fasta_path,
                target_file=str(bed_path),
                itd_length_cutoff=10,
                allowed_mismatches_for_sr_rescue=1,
                allowed_mismatches_for_insertion=2,
                logger=logger,
                max_depth=-1,
            )
