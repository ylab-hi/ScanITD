# Output Format

ScanITD writes results in **VCF 4.3** format.

## File structure

A ScanITD VCF file contains:

1. **Meta-information lines** (starting with `##`)
2. **Header line** (starting with `#CHROM`)
3. **Data records** — one per detected event

## VCF columns

| Column | Field | Description |
|--------|-------|-------------|
| 1 | `CHROM` | Chromosome of the event |
| 2 | `POS` | 1-based reference position |
| 3 | `ID` | Sequential event identifier |
| 4 | `REF` | Reference base at `POS` |
| 5 | `ALT` | `<TDUP>` or `<INS>` symbolic allele |
| 6 | `QUAL` | `.` (not computed) |
| 7 | `FILTER` | `.` (no filter applied) |
| 8 | `INFO` | Semicolon-separated key=value fields (see below) |
| 9 | `FORMAT` | `GT` |
| 10 | `SAMPLE` | Genotype call (`0/1`) |

## INFO fields

| Key | Type | Description |
|-----|------|-------------|
| `SVTYPE` | String | Event type: `TDUP` (tandem duplication) or `INS` (insertion) |
| `SVLEN` | Integer | Length of the duplicated / inserted sequence (bp) |
| `END` | Integer | 1-based end coordinate of the event |
| `CHR2` | String | Chromosome of the end coordinate (always same as CHROM) |
| `OAO` | Integer | **Original** alternate allele observations (SA-tag reads only) |
| `AO` | Integer | **Rescued** alternate allele observations (SA-tag + rescued soft-clips) |
| `DP` | Integer | Total read depth at the locus |
| `AF` | Float | Variant allele frequency (`AO / DP`) |
| `SEQ` | String | Duplicated / inserted sequence |
| `INSSEQ` | String | Microinsertion sequence at the breakpoint (`.` if none) |
| `HOMSEQ` | String | Microhomology sequence at the breakpoint (`.` if none) |
| `SVMETHOD` | String | Detection method (`ScanITD2`) |

## Event types

### TDUP — Tandem Duplication (ITD)

An internal tandem duplication where a segment of the genome is duplicated
in tandem at the same locus.

```
REF:  ─────[ABCDE]──────────────────
ALT:  ─────[ABCDE][ABCDE]───────────
                   ↑ duplicated copy
POS = start of the duplicated segment
END = POS + SVLEN
```

### INS — Novel Insertion

A large insertion in the CIGAR string that does **not** match a self-loop
(tandem duplication) pattern.

```
REF:  ─────────────────────
ALT:  ──────[NOVEL SEQ]────
POS = insertion site
END = POS (zero-length on reference)
```

## Example records

```
##fileformat=VCFv4.3
##fileDate=20241108
##source=ScanITDv0.9.1
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The type of event, TDUP, INS.">
##INFO=<ID=AO,Number=1,Type=Integer,Description="Alternate allele observations">
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency ...">
...
#CHROM  POS       ID  REF  ALT     QUAL  FILTER  INFO                                          FORMAT  sample
chr13   28034008  1   A    <TDUP>  .     .       SVTYPE=TDUP;OAO=45;AO=58;DP=120;AF=0.483;... GT      0/1
chr13   28033998  2   G    <INS>   .     .       SVTYPE=INS;OAO=6;AO=6;DP=110;AF=0.0545;...   GT      0/1
```

## Breakpoint region fields

The `INSSEQ` and `HOMSEQ` fields describe what lies at the ITD junction:

| Scenario | `INSSEQ` | `HOMSEQ` |
|----------|----------|----------|
| Blunt-end junction | `.` | `.` |
| Microinsertion (novel bases between copies) | inserted sequence | `.` |
| Microhomology (overlapping bases at junction) | `.` | overlapping sequence |
