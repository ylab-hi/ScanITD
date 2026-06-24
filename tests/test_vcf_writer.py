"""Tests for scanitd.writer.vcf_writer — VCFWriter, feature helpers."""

import pytest

from scanitd.base import Event, MicroRegion
from scanitd.writer.vcf_writer import (
    VCFWriter,
    get_vcf_features_from_event,
    obtain_reference_from_bam_header,
    vcf_feature_transformer,
)


# ---------------------------------------------------------------------------
# Minimal BAM header fixture (module-level, not from conftest to keep
# writer tests self-contained)
# ---------------------------------------------------------------------------

MINIMAL_BAM_HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [
        {"SN": "chr1", "LN": 248956422},
        {"SN": "chr2", "LN": 242193529},
    ],
    "PG": [{"ID": "bwa", "CL": "bwa mem ref.fa sample.fastq"}],
}


def make_event(event_type="TDUP", chrom="chr1", ref_start=100, size=50,
               seq="ACGT", oao=5, ao=8, dp=40,
               ref_allele="A", alt_allele="<TDUP>",
               micro=""):
    region = MicroRegion(micro)
    return Event.new(
        event_type=event_type,
        event_id=(chrom, ref_start, size, seq, region),
        oao=oao,
        ao=ao,
        dp=dp,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
    )


# ---------------------------------------------------------------------------
# obtain_reference_from_bam_header
# ---------------------------------------------------------------------------

class TestObtainReferenceFromBamHeader:
    def test_known_aligner_returns_command_line(self):
        result = obtain_reference_from_bam_header(MINIMAL_BAM_HEADER)
        assert result == "bwa mem ref.fa sample.fastq"

    def test_unknown_aligner_returns_unknown(self):
        header = {"PG": [{"ID": "custom_aligner", "CL": "custom ..."}]}
        assert obtain_reference_from_bam_header(header) == "Unknown"

    def test_no_pg_section_returns_unknown(self):
        assert obtain_reference_from_bam_header({}) == "Unknown"

    def test_case_insensitive_aligner_match(self):
        header = {"PG": [{"ID": "BWA", "CL": "BWA mem ..."}]}
        assert obtain_reference_from_bam_header(header) == "BWA mem ..."


# ---------------------------------------------------------------------------
# get_vcf_features_from_event
# ---------------------------------------------------------------------------

class TestGetVcfFeaturesFromEvent:
    def test_tdup_features(self):
        ev = make_event("TDUP")
        feat = get_vcf_features_from_event(ev)
        assert feat["SVTYPE"] == "TDUP"
        assert feat["CHROM"] == "chr1"
        # POS is 1-based
        assert feat["POS"] == "101"
        # END is 1-based (ref_start + event_size + 1)
        assert feat["END"] == str(ev.end + 1)
        assert feat["SEQ"] == "ACGT"
        assert feat["INSSEQ"] == "."   # blunt end
        assert feat["HOMSEQ"] == "."

    def test_ins_features(self):
        ev = make_event("INS", seq="GCGCGC", ref_start=200, size=6)
        feat = get_vcf_features_from_event(ev)
        assert feat["SVTYPE"] == "INS"
        assert feat["SEQ"] == "GCGCGC"

    def test_microinsertion_populates_insseq(self):
        ev = make_event(micro="+TTT")
        feat = get_vcf_features_from_event(ev)
        assert feat["INSSEQ"] == "TTT"
        assert feat["HOMSEQ"] == "."

    def test_microhomology_populates_homseq(self):
        ev = make_event(micro="-GGG")
        feat = get_vcf_features_from_event(ev)
        assert feat["HOMSEQ"] == "GGG"
        assert feat["INSSEQ"] == "."

    def test_af_formatted_as_3_sig_figs(self):
        ev = make_event(ao=1, dp=3)
        feat = get_vcf_features_from_event(ev)
        # 1/3 ≈ 0.333 formatted with .3g
        assert "0.333" in feat["AF"] or "0.33" in feat["AF"]


# ---------------------------------------------------------------------------
# vcf_feature_transformer
# ---------------------------------------------------------------------------

class TestVcfFeatureTransformer:
    def _make_feature_dict(self, **overrides):
        base = {
            "CHROM": "chr1", "POS": "101", "REF": "A",
            "ALT": "<TDUP>", "SVTYPE": "TDUP",
            "CHR2": "chr1", "END": "151",
            "OAO": "5", "AO": "8", "DP": "40", "AF": "0.2",
            "SVLEN": "50", "SVMETHOD": "ScanITD2",
            "SEQ": "ACGT", "INSSEQ": ".", "HOMSEQ": ".",
        }
        base.update(overrides)
        return base

    def test_returns_10_fields(self):
        row = vcf_feature_transformer(self._make_feature_dict(), 1)
        assert len(row) == 10

    def test_field_order(self):
        row = vcf_feature_transformer(self._make_feature_dict(), 42)
        chrom, pos, id_, ref, alt, qual, filt, info, fmt, sample = row
        assert chrom == "chr1"
        assert pos == "101"
        assert id_ == "42"
        assert qual == "."
        assert filt == "."
        assert fmt == "GT"
        assert sample == "0/1"

    def test_info_field_contains_svtype(self):
        row = vcf_feature_transformer(self._make_feature_dict(), 1)
        info = row[7]
        assert "SVTYPE=TDUP" in info

    def test_event_id_as_string(self):
        """event_id may be passed as a str (from write_events_to_vcf)."""
        row = vcf_feature_transformer(self._make_feature_dict(), "7")
        assert row[2] == "7"


# ---------------------------------------------------------------------------
# VCFWriter  (integration — writes a real file)
# ---------------------------------------------------------------------------

class TestVCFWriter:
    def test_write_header_on_open(self, tmp_path):
        out = tmp_path / "out.vcf"
        writer = VCFWriter(str(out), MINIMAL_BAM_HEADER)
        with writer.open():
            pass
        content = out.read_text()
        assert content.startswith("##fileformat=VCFv4.3")
        assert "##contig=<ID=chr1" in content
        assert "#CHROM\tPOS\tID\tREF\tALT" in content

    def test_sample_name_from_stem(self, tmp_path):
        out = tmp_path / "mysample.vcf"
        writer = VCFWriter(str(out), MINIMAL_BAM_HEADER)
        assert writer.sample_name == "mysample"

    # Bug #9 regression -------------------------------------------------------
    def test_write_data_actually_writes_event(self, tmp_path):
        """write_data(Event, id) must produce a VCF data line — not be a no-op."""
        out = tmp_path / "test.vcf"
        ev = make_event("TDUP")
        writer = VCFWriter(str(out), MINIMAL_BAM_HEADER)
        with writer.open():
            writer.write_data(ev, "1")
        content = out.read_text()
        lines = [l for l in content.splitlines() if not l.startswith("#")]
        assert len(lines) == 1, "Expected exactly one data line in the VCF"
        fields = lines[0].split("\t")
        assert fields[0] == "chr1"
        assert "SVTYPE=TDUP" in fields[7]

    def test_write_multiple_events(self, tmp_path):
        out = tmp_path / "multi.vcf"
        events = [make_event("TDUP", ref_start=i * 100) for i in range(3)]
        writer = VCFWriter(str(out), MINIMAL_BAM_HEADER)
        with writer.open():
            for idx, ev in enumerate(events, 1):
                writer.write_data(ev, str(idx))
        lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        assert len(lines) == 3

    def test_write_data_unknown_type_raises(self, tmp_path):
        """write_data with an unsupported type must raise NotImplementedError."""
        out = tmp_path / "err.vcf"
        writer = VCFWriter(str(out), MINIMAL_BAM_HEADER)
        with writer.open():
            with pytest.raises(NotImplementedError):
                writer.write_data("not_an_event", "1")

    def test_is_opened_lifecycle(self, tmp_path):
        out = tmp_path / "lifecycle.vcf"
        writer = VCFWriter(str(out), MINIMAL_BAM_HEADER)
        assert not writer.is_opened
        with writer.open():
            assert writer.is_opened
        assert not writer.is_opened
