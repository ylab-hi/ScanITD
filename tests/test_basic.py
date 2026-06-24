"""Tests for scanitd.base.basic — Event, Interval, Intervals, MicroRegion, Strand, MappingMode."""

import pytest

from scanitd.base import Event, Interval, Intervals, MappingMode, MicroRegion, Strand


# ---------------------------------------------------------------------------
# MicroRegion
# ---------------------------------------------------------------------------

class TestMicroRegion:
    def test_blunt_end(self):
        m = MicroRegion("")
        assert m.micro_type == "blunt_end"
        assert m.sequence == ""
        assert m.length == 0

    def test_microinsertion(self):
        m = MicroRegion("+ACGT")
        assert m.micro_type == "microinsertion"
        assert m.sequence == "ACGT"
        assert m.length == 4

    def test_microhomology(self):
        m = MicroRegion("-TTG")
        assert m.micro_type == "microhomology"
        assert m.sequence == "TTG"
        assert m.length == 3

    def test_equality(self):
        assert MicroRegion("+AC") == MicroRegion("+AC")
        assert MicroRegion("+AC") != MicroRegion("+AG")

    def test_hash_consistent_with_eq(self):
        a, b = MicroRegion("-GC"), MicroRegion("-GC")
        assert hash(a) == hash(b)


# ---------------------------------------------------------------------------
# Strand
# ---------------------------------------------------------------------------

class TestStrand:
    def test_from_str_forward(self):
        assert Strand.from_str("+") == Strand.Forward

    def test_from_str_reverse(self):
        assert Strand.from_str("-") == Strand.Reverse

    def test_from_str_passthrough(self):
        """Passing an existing Strand instance returns it unchanged."""
        s = Strand.Forward
        assert Strand.from_str(s) is s

    def test_invalid_strand(self):
        with pytest.raises(ValueError, match="Invalid strand"):
            Strand.from_str("x")

    def test_is_forward_reverse(self):
        assert Strand.Forward.is_forward()
        assert not Strand.Forward.is_reverse()
        assert Strand.Reverse.is_reverse()

    def test_reversed(self):
        assert Strand.Forward.reversed() == Strand.Reverse
        assert Strand.Reverse.reversed() == Strand.Forward

    def test_str(self):
        assert str(Strand.Forward) == "+"
        assert str(Strand.Reverse) == "-"


# ---------------------------------------------------------------------------
# MappingMode
# ---------------------------------------------------------------------------

class TestMappingMode:
    def test_from_int(self):
        assert MappingMode.from_int(1) == MappingMode.MS
        assert MappingMode.from_int(2) == MappingMode.SM

    def test_from_int_passthrough(self):
        assert MappingMode.from_int(MappingMode.SM) == MappingMode.SM

    def test_invalid_int(self):
        with pytest.raises(ValueError, match="Invalid mode"):
            MappingMode.from_int(99)

    def test_reversed(self):
        assert MappingMode.SM.reversed() == MappingMode.MS
        assert MappingMode.MS.reversed() == MappingMode.SM

    def test_predicates(self):
        assert MappingMode.SM.is_sm()
        assert not MappingMode.SM.is_ms()
        assert MappingMode.MS.is_ms()


# ---------------------------------------------------------------------------
# Event
# ---------------------------------------------------------------------------

class TestEvent:
    def test_tdup_event_end(self, tdup_event):
        """TDUP end = ref_start + event_size."""
        assert tdup_event.end == tdup_event.ref_start + tdup_event.event_size

    def test_ins_event_end(self, ins_event):
        """INS end = ref_start (zero-length on reference)."""
        assert ins_event.end == ins_event.ref_start

    def test_af_computed(self, tdup_event):
        expected = round(tdup_event.ao / tdup_event.dp, 4)
        assert tdup_event.af == expected

    def test_chrom2_equals_chrom(self, tdup_event):
        assert tdup_event.chrom2 == tdup_event.chrom

    # Bug #2 regression -------------------------------------------------------
    def test_unknown_event_type_raises(self, blunt_region):
        """Event.new() with an unknown type must raise ValueError, not UnboundLocalError."""
        with pytest.raises(ValueError, match="Unknown event_type"):
            Event.new(
                event_type="UNKNOWN",
                event_id=("chr1", 0, 10, "ACGT", blunt_region),
                oao=1,
                ao=1,
                dp=10,
            )

    def test_hash_and_repr_do_not_crash(self, tdup_event):
        assert isinstance(hash(tdup_event), int)
        assert "Event(" in repr(tdup_event)


# ---------------------------------------------------------------------------
# Interval
# ---------------------------------------------------------------------------

class TestInterval:
    def test_basic_construction(self, simple_interval):
        assert simple_interval.start == 10
        assert simple_interval.end == 20

    def test_len(self, simple_interval):
        assert len(simple_interval) == 10

    def test_start_gt_end_raises(self):
        with pytest.raises(ValueError):
            Interval(20, 10)

    def test_contains(self):
        outer = Interval(0, 100)
        inner = Interval(10, 50)
        assert inner in outer
        assert not outer in inner

    def test_getitem(self, simple_interval):
        assert simple_interval[0] == 10
        assert simple_interval[1] == 20

    def test_getitem_out_of_range(self, simple_interval):
        with pytest.raises((ValueError, IndexError)):
            _ = simple_interval[2]

    # Bug #3 regression -------------------------------------------------------
    def test_setitem_valid_does_not_raise(self, simple_interval):
        """__setitem__ must NOT raise for valid indices 0 and 1."""
        simple_interval[0] = 5
        assert simple_interval.start == 5
        simple_interval[1] = 25
        assert simple_interval.end == 25

    def test_setitem_invalid_raises(self, simple_interval):
        with pytest.raises(IndexError):
            simple_interval[2] = 99

    def test_add_int(self, simple_interval):
        shifted = simple_interval + 5
        assert shifted.start == 15
        assert shifted.end == 25

    def test_sub_int(self, simple_interval):
        shifted = simple_interval - 5
        assert shifted.start == 5
        assert shifted.end == 15

    def test_equality(self):
        assert Interval(1, 5) == Interval(1, 5)
        assert Interval(1, 5) != Interval(1, 6)

    def test_hash_consistent(self):
        assert hash(Interval(1, 5)) == hash(Interval(1, 5))

    def test_iter(self, simple_interval):
        start, end = simple_interval
        assert start == 10
        assert end == 20

    def test_from_list(self):
        iv = Interval.from_list([3, 7])
        assert iv.start == 3 and iv.end == 7

    def test_repr(self, simple_interval):
        assert "[10, 20)" in repr(simple_interval)


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------

class TestIntervals:
    def test_getitem(self, two_exon_intervals):
        assert two_exon_intervals[0] == Interval(0, 10)
        assert two_exon_intervals[1] == Interval(20, 30)

    def test_len_returns_count_not_span(self, two_exon_intervals):
        """len() returns count of exons, not genomic span."""
        assert len(two_exon_intervals) == 2

    def test_contains(self, two_exon_intervals):
        assert Interval(0, 10) in two_exon_intervals
        assert Interval(0, 11) not in two_exon_intervals

    def test_iter(self, two_exon_intervals):
        items = list(two_exon_intervals)
        assert items == [Interval(0, 10), Interval(20, 30)]

    # Bug #4 regression -------------------------------------------------------
    def test_reverse_mutates_list(self, two_exon_intervals):
        """reverse() must actually mutate self.exon_list in place."""
        original_first = two_exon_intervals[0]
        two_exon_intervals.reverse()
        assert two_exon_intervals[0] != original_first
        assert two_exon_intervals[0] == Interval(20, 30)

    def test_reversed_returns_new_object(self, two_exon_intervals):
        rev = two_exon_intervals.reversed()
        assert rev[0] == Interval(20, 30)
        # original unchanged
        assert two_exon_intervals[0] == Interval(0, 10)

    def test_add_int(self, two_exon_intervals):
        shifted = two_exon_intervals + 5
        assert shifted[0] == Interval(5, 15)
        assert shifted[1] == Interval(25, 35)

    def test_append_interval(self, two_exon_intervals):
        two_exon_intervals.append(Interval(40, 50))
        assert len(two_exon_intervals) == 3

    def test_append_tuple(self, two_exon_intervals):
        two_exon_intervals.append((40, 50))
        assert two_exon_intervals[-1] == Interval(40, 50)

    def test_introns(self, two_exon_intervals):
        introns = two_exon_intervals.introns()
        assert introns is not None
        assert len(introns) == 1
        assert introns[0] == Interval(10, 20)

    def test_introns_single_exon_returns_none(self):
        ivs = Intervals([Interval(0, 10)])
        assert ivs.introns() is None

    def test_first_last(self, two_exon_intervals):
        assert two_exon_intervals.first == Interval(0, 10)
        assert two_exon_intervals.last == Interval(20, 30)

    def test_from_list(self):
        ivs = Intervals.from_list([[0, 5], [10, 15]])
        assert ivs[0] == Interval(0, 5)
