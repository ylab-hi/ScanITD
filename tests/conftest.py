"""Shared pytest fixtures for ScanITD tests."""

import pytest

from scanitd.base import Event, Interval, Intervals, MicroRegion


# ---------------------------------------------------------------------------
# MicroRegion fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def blunt_region():
    return MicroRegion("")


@pytest.fixture
def insertion_region():
    return MicroRegion("+ACGT")


@pytest.fixture
def homology_region():
    return MicroRegion("-TTG")


# ---------------------------------------------------------------------------
# Event fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tdup_event(blunt_region):
    return Event.new(
        event_type="TDUP",
        event_id=("chr1", 100, 50, "ACGTACGT", blunt_region),
        oao=5,
        ao=8,
        dp=40,
        ref_allele="A",
        alt_allele="<TDUP>",
    )


@pytest.fixture
def ins_event(blunt_region):
    return Event.new(
        event_type="INS",
        event_id=("chr2", 200, 12, "GCTAGCTAGCTA", blunt_region),
        oao=3,
        ao=4,
        dp=20,
        ref_allele="G",
        alt_allele="GGCTAGCTAGCTA",
    )


# ---------------------------------------------------------------------------
# Interval / Intervals fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_interval():
    return Interval(10, 20)


@pytest.fixture
def two_exon_intervals():
    return Intervals([Interval(0, 10), Interval(20, 30)])
