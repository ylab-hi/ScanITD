"""Fundamental data types and enumerations for ScanITD.

Provides genomic interval arithmetic, structural variant event containers,
alignment mode enums, CIGAR operation codes, and strand representation.
"""

from __future__ import annotations

from enum import Enum, IntEnum


class MicroRegion:
    """Store microinsertion or microhomology event."""

    __slots__ = (
        "length",
        "micro_type",
        "sequence",
    )

    def __init__(
        self,
        input_sequence: str,
    ) -> None:
        """Initialize a MicroRegion from a prefixed sequence string.

        Args:
            input_sequence: Breakpoint region sequence with a prefix indicating type:
                ``+<seq>`` for microinsertion, ``-<seq>`` for microhomology,
                or empty string / plain sequence for blunt-end.
        """
        if input_sequence.startswith("+"):
            self.micro_type = "microinsertion"
            self.sequence = input_sequence[1:]
        elif input_sequence.startswith("-"):
            self.micro_type = "microhomology"
            self.sequence = input_sequence[1:]
        else:
            self.micro_type = "blunt_end"
            self.sequence = ""
        self.length = len(self.sequence)

    def __hash__(self) -> int:
        """Get the hash value of the event.
        :return: hash value of the event
        """
        return hash(self.micro_type) ^ hash(self.sequence) ^ hash(self.length)

    def __eq__(self, other) -> bool:
        """Return True if both MicroRegion objects share the same type, sequence and length."""
        if not isinstance(other, MicroRegion):
            return False

        # Compare `name` and `age` attributes for equality
        return self.micro_type == other.micro_type and self.sequence == other.sequence and self.length == other.length

    def __repr__(self) -> str:
        """Get the representation of the event.

        :return: representation of the event
        """
        return f"MicroRegion({self.micro_type=}, {self.sequence=} {self.length=})"


class Event:
    """Store TDUP or INS event."""

    __slots__ = (
        "af",
        "alt_allele",
        "oao",
        "ao",
        "break_point_region",
        "chrom",
        "chrom2",
        "dp",
        "end",
        "event_sequence",
        "event_size",
        "event_type",
        "ref_allele",
        "ref_start",
    )

    def __init__(
        self,
        chrom: str,
        ref_start: int,
        event_size: int,
        event_sequence: str,
        event_type: str,
        oao: int,
        ao: int,
        dp: int,
        end: int,
        ref_allele: str | None = None,
        alt_allele: str | None = None,
        break_point_region: MicroRegion | None = None,
    ) -> None:
        """Initialize an Event.

        Args:
            chrom: Chromosome name (e.g. ``chr1``).
            ref_start: 0-based start position on the reference.
            event_size: Size of the structural variant in base pairs.
            event_sequence: Sequence of the event (inserted or duplicated bases).
            event_type: Event class; one of ``TDUP`` or ``INS``.
            oao: Original (pre-rescue) alternate allele observation count.
            ao: Rescued alternate allele observation count.
            dp: Total read depth at the locus.
            end: 0-based end position on the reference.
            ref_allele: Reference base at the locus, or ``None``.
            alt_allele: Alternate allele sequence, or ``None``.
            break_point_region: Micro-homology / micro-insertion at the breakpoint,
                or ``None`` for blunt-end junctions.
        """
        self.chrom = chrom
        self.ref_start = ref_start
        self.event_size = event_size
        self.event_sequence = event_sequence
        self.event_type = event_type
        self.oao = oao
        self.ao = ao
        self.dp = dp
        self.af = round(float(ao / dp), 4)
        self.chrom2 = self.chrom
        self.end = end
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.break_point_region = break_point_region

    def __hash__(self) -> int:
        """Get the hash value of the event.

        :return: hash value of the event
        """
        return hash(self.chrom) ^ hash(self.ref_start) ^ hash(self.event_size) ^ hash(self.event_sequence) ^ hash(self.event_type)

    def __repr__(self) -> str:
        """Get the representation of the event.

        :return: representation of the event
        """
        return f"Event({self.chrom=}, {self.ref_start=}, {self.event_size=}, {self.event_sequence=}, {self.event_type=} " f"{self.ref_allele=}, {self.alt_allele=})"

    @classmethod
    def new(
        cls,
        event_type: str,
        event_id: tuple[str, int, int, str, MicroRegion],
        oao: int,
        ao: int,
        dp: int,
        ref_allele: str | None = None,
        alt_allele: str | None = None,
    ):
        """Construct an Event from a structured event ID tuple.

        Args:
            event_type: Either ``TDUP`` or ``INS``.
            event_id: A 5-tuple of (chrom, ref_start, event_size, event_sequence,
                break_point_region).
            oao: Original (pre-rescue) alternate allele observation count.
            ao: Rescued alternate allele observation count.
            dp: Total read depth at the locus.
            ref_allele: Reference base at the locus.
            alt_allele: Alternate allele sequence.

        Returns:
            A new Event instance.

        Raises:
            ValueError: If event_type is not ``TDUP`` or ``INS``.
        """
        chrom, ref_start, event_size, event_sequence, break_point_region = event_id

        if event_type == "TDUP":
            end = ref_start + event_size
        elif event_type == "INS":
            end = ref_start
        else:
            msg = f"Unknown event_type: {event_type!r}"
            raise ValueError(msg)

        return cls(
            chrom,
            ref_start,
            event_size,
            event_sequence,
            event_type,
            oao,
            ao,
            dp,
            end,
            ref_allele,
            alt_allele,
            break_point_region,
        )


class MappingMode(IntEnum):
    """Chimeric-read alignment mode relative to the soft-clipped end.

    - MS (1): Match then Softclip — the matched segment precedes the softclip.
    - SM (2): Softclip then Match — the softclip precedes the matched segment.
    """

    Type0 = 0
    MS = 1  # 1 MS
    SM = 2  # 2 SM

    # fmt: off
    def reversed(self):
        """Return the opposite mapping mode (MS <-> SM)."""
        return self.SM if self == self.MS else self.MS

    def is_sm(self):
        """Return True if this mode is SM (Softclip-Match)."""
        return self == self.SM

    def is_ms(self):
        """Return True if this mode is MS (Match-Softclip)."""
        return self == self.MS

    def to_str(self):
        """Return the string representation (``SM`` or ``MS``)."""
        return "SM" if self == self.SM else "MS"
    # fmt: on

    @classmethod
    def from_int(cls, mode: int | MappingMode):
        """Construct a MappingMode from an integer code.

        Args:
            mode: 0 for Type0, 1 for MS, 2 for SM, or an existing MappingMode.

        Returns:
            Corresponding MappingMode member.

        Raises:
            ValueError: If mode is not 0, 1, or 2.
        """
        if isinstance(mode, cls):
            return mode

        if mode == 0:
            return cls.Type0
        if mode == 1:
            return cls.MS
        if mode == 2:  # noqa: PLR2004
            return cls.SM

        msg = f"Invalid mode: {mode}"
        raise ValueError(msg)


class AnnotationCode(IntEnum):
    """Integer codes used to annotate chimeric-read breakpoint types."""

    Type1 = 1
    Type2 = 2

    def reversed(self):
        """Return the other AnnotationCode (Type1 <-> Type2)."""
        return self.Type2 if self == self.Type1 else self.Type1


class Strand(Enum):
    """Strand."""

    Forward = "+"
    Reverse = "-"

    # fmt: off
    def is_reverse(self):
        """Return True if this strand is the reverse strand."""
        return self == Strand.Reverse

    def is_forward(self):
        """Return True if this strand is the forward strand."""
        return self == Strand.Forward

    def reversed(self) -> Strand:
        """Return the opposite strand."""
        return Strand.Reverse if self.is_forward() else Strand.Forward
    # fmt: on

    def __str__(self) -> str:
        if self.is_forward():
            return "+"
        return "-"

    @classmethod
    def from_str(cls, strand: str | Strand):
        """Construct a Strand from a string or pass through an existing Strand.

        Args:
            strand: ``+`` for forward, ``-`` for reverse, or an existing Strand.

        Returns:
            Corresponding Strand member.

        Raises:
            ValueError: If strand string is not ``+`` or ``-``.
        """
        if isinstance(strand, cls):
            return strand
        if strand == "+":
            return cls.Forward
        if strand == "-":
            return cls.Reverse
        msg = f"Invalid strand: {strand}"
        raise ValueError(msg)


# // #define BAM_CMATCH      0
# // #define BAM_CINS        1
# // #define BAM_CDEL        2
# // #define BAM_CREF_SKIP   3
# // #define BAM_CSOFT_CLIP  4
# // #define BAM_CHARD_CLIP  5
# // #define BAM_CEQUAL      7
# // #define BAM_CPAD        6
# // #define BAM_CDIFF       8
# // #define BAM_CBACK       9
# // #define BAM_CIGAR_STR   "MIDNSHP=XB"


class CigarCode(IntEnum):
    """Integer enum mapping CIGAR operation characters to their SAM/BAM encoding codes as defined in the SAM specification."""

    Match = 0
    Insertion = 1
    Del = 2
    Ref_skip = 3
    Soft_clip = 4
    Hard_clip = 5
    Pad = 6
    Equal = 7
    Diff = 8
    Back = 9


class Interval:
    """0-based interval is used to represent the interval of genomics.

    .. note:
        start <= end

    """

    start: int
    end: int
    _index: int = 2

    def __init__(self, start: int, end: int):
        """Initialize Interval."""
        if start > end:
            msg = f"start: {start} > end: {end}"
            raise ValueError(msg)

        self.start = start
        self.end = end

    def __eq__(self, other: Interval):
        """Return True if both intervals share the same start and end positions."""
        if isinstance(other, Interval):
            return self.start == other.start and self.end == other.end
        return False

    def __contains__(self, item: Interval):
        """Check if item is in the interval."""
        return self.start <= item.start and item.end <= self.end

    def __hash__(self) -> int:
        """Return hash based on start and end."""
        return hash(self.start) ^ hash(self.end)

    def __add__(self, other: Interval | int) -> Interval:
        """Return a new Interval shifted by an integer offset or element-wise sum with another Interval."""
        if isinstance(other, int):
            return Interval(self.start + other, self.end + other)

        if isinstance(other, Interval):
            return Interval(self.start + other.start, self.end + other.end)

        msg = f"{other} is not Interval or int"
        raise TypeError(msg)

    def __sub__(self, other: Interval | int) -> Interval:
        """Return a new Interval shifted by an integer offset or element-wise difference with another Interval."""
        if isinstance(other, int):
            return Interval(self.start - other, self.end - other)

        if isinstance(other, Interval):
            return Interval(self.start - other.start, self.end - other.end)

        msg = f"{other} is not Interval or int"
        raise TypeError(msg)

    def __setitem__(self, index: int, value: int):
        """Set start (index 0) or end (index 1)."""
        if index == 0:
            self.start = value
        elif index == 1:
            self.end = value
        else:
            msg = f"index: {index} is out"
            raise IndexError(msg)

    def __getitem__(self, index: int):
        """Get item from interval."""
        if index >= Interval._index:
            msg = f"index: {index} is out"
            raise ValueError(msg)

        if index == 0:
            return self.start
        if index == 1:
            return self.end
        msg = f"index: {index} is out"
        raise IndexError(msg)

    def __len__(self):
        """Return the length of the interval (end - start)."""
        return self.end - self.start

    def __repr__(self) -> str:
        return f"[{self.start}, {self.end})"

    def __iter__(self):
        """Iterate over (start, end)."""
        return iter((self.start, self.end))

    __str__ = __repr__

    @classmethod
    def from_list(cls, item: list[int] | tuple[int, int]):
        """Create Interval from tuple."""
        return cls(*item)

    def join(self, other: Interval) -> tuple[Interval, Interval] | None:
        """Compute the intersection and union of two overlapping intervals.

        Args:
            other: Another Interval to join with.

        Returns:
            A tuple of (overlap, union) Intervals, or None if the intervals do not overlap.

        Raises:
            ValueError: If other is not an Interval.
        """
        if isinstance(other, Interval):
            if self.start <= other.start < self.end:
                if other.end < self.end:
                    # s |o o| s
                    return Interval(other.start, other.end), Interval(
                        self.start,
                        self.end,
                    )
                # s |o s| o
                return Interval(other.start, self.end), Interval(self.start, other.end)

            if other.start <= self.start < other.end:
                # o |s s| o
                if self.end < other.end:
                    return Interval(self.start, self.end), Interval(
                        other.start,
                        other.end,
                    )
                # o |s o| s
                return Interval(self.start, other.end), Interval(other.start, self.end)

            if self.start == other.end:
                #  o  os  s
                return Interval(other.end, self.start), Interval(other.start, self.end)

            if self.end == other.start:
                # s  so  o
                return Interval(self.end, other.start), Interval(
                    self.start,
                    other.end,
                )

            # no overlap o  o  s  s or s  s  o  o
            return None

        msg = f"{other} is not a Interval"
        raise ValueError(msg)

    def contain(self, other: Interval, *, same_left=False, same_right=False) -> bool:
        """Check whether this interval contains another, with optional boundary matching.

        Args:
            other: The Interval to test for containment.
            same_left: If True, left boundaries must match exactly.
            same_right: If True, right boundaries must match exactly.

        Returns:
            True if other is contained within self under the specified boundary conditions.

        Raises:
            ValueError: If other is not an Interval.
        """
        if isinstance(other, Interval):
            if not same_left and not same_right:
                return other in self
            if same_left and same_right:
                return other == self
            if same_left:
                return self.end >= other.end
            if same_right:
                return self.start <= other.start

        msg = f"{other} is not Interval"
        raise ValueError(msg)


class Intervals:
    """A sequence of Interval objects representing non-overlapping genomic regions (e.g. exons).

    Args:
        exon_list: Ordered list of Interval objects.
    """

    def __init__(self, exon_list: list[Interval]) -> None:
        self.exon_list = exon_list

    # fmt: off
    def __hash__(self) -> int: return hash(tuple(self.exon_list))
    def __getitem__(self, index: int) -> Interval: return self.exon_list[index]
    def __setitem__(self, index: int, value: Interval) -> None: self.exon_list[index] = value
    def __len__(self): return len(self.exon_list)
    def __contains__(self, item: Interval): return any(item == exon for exon in self.exon_list)
    def __iter__(self): return iter(self.exon_list)
    def __repr__(self) -> str: return f"Exons({self.exon_list})"
    def __str__(self) -> str: return "[" + ",".join([f"{interval.start}-{interval.end}" for interval in self.exon_list]) + "]"
    def reverse(self) -> None:
        """Reverse the order of intervals in-place."""
        self.exon_list.reverse()

    def reversed(self) -> Intervals:
        """Return a new Intervals with the order reversed."""
        return Intervals(self.exon_list[::-1])

    @property
    def first(self) -> Interval:
        """Return the first interval."""
        return self.exon_list[0]

    @property
    def last(self) -> Interval:
        """Return the last interval."""
        return self.exon_list[-1]
    # fmt: on

    def __add__(self, other: int | Intervals) -> Intervals:
        """Shift all intervals by an integer offset, or add element-wise with another Intervals."""
        if isinstance(other, int):
            return Intervals([exon + other for exon in self.exon_list])

        if isinstance(other, Intervals):
            return Intervals(
                [exon + other_exon for exon, other_exon in zip(self.exon_list, other.exon_list, strict=False)],
            )

        message = f"{other} is not int or Intervals"
        raise TypeError(message)

    def __sub__(self, other: int | Intervals) -> Intervals:
        """Shift all intervals by an integer offset, or subtract element-wise with another Intervals."""
        if isinstance(other, int):
            return Intervals([exon - other for exon in self.exon_list])

        if isinstance(other, Intervals):
            return Intervals(
                [exon - other_exon for exon, other_exon in zip(self.exon_list, other.exon_list, strict=False)],
            )

        message = f"{other} is not int or Intervals"
        raise TypeError(message)

    def __eq__(self, other: Intervals) -> bool:
        """Return True if all constituent intervals are equal."""
        if isinstance(other, Intervals):
            return self.exon_list == other.exon_list
        return False

    def append(self, item: Interval | tuple[int, int]) -> None:
        """Append item to exon_list."""
        if isinstance(item, tuple):
            self.exon_list.append(Interval(*item))
        elif isinstance(item, Interval):
            self.exon_list.append(item)
        else:
            msg = f"item: {item} is not Interval or tuple"
            raise TypeError(msg)

    def extend(self, item: Intervals) -> None:
        """Extend exon_list with item."""
        if isinstance(item, Intervals):
            self.exon_list.extend(item.exon_list)
        else:
            msg = f"item: {item} is not Intervals"
            raise TypeError(msg)

    @classmethod
    def from_list(cls, item: list[list[int] | tuple[int, int]]):
        """Create Exons from list."""
        return cls(exon_list=[Interval.from_list(exon) for exon in item])

    def introns(self) -> Intervals | None:
        """Compute the intron intervals between consecutive exons.

        Returns:
            An Intervals containing the gaps between exons, or None if fewer than
            2 exons exist.
        """
        if len(self) < Interval._index:
            return None

        introns = Intervals([])
        for exon_group in zip(self.exon_list, self.exon_list[1:], strict=False):
            introns.append(Interval(exon_group[0].end, exon_group[1].start))

        return introns


Exon = Interval
Exons = Intervals
Introns = Intervals
