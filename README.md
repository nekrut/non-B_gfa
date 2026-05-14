# non-B-gfa

Tools for finding non-B DNA forming motifs in nucleic acid sequences.
Developed at NCI-Frederick / Frederick National Laboratory.

DNA exists in many possible conformations beyond the canonical right-handed
B-form Watson-Crick double helix: cruciforms, triplex, slipped (hairpin)
structures, G-quadruplexes, left-handed Z-DNA, and others. Several
publications have provided significant evidence that non-B DNA structures
may play a role in DNA instability and mutagenesis, leading to both DNA
rearrangements and increased mutational rates — hallmarks of cancer.

**Web submission:** https://nonb-abcc.ncifcrf.gov/apps/site/default

The website results use the default values of `gfa` and should match the
example outputs included in `test_files.tar` (see Quick Start).

**Please cite:** Non-B DB v2.0: a database of predicted non-B DNA-forming
motifs and its associated tools. Regina Z. Cer et al.
*Nucl. Acids Res.* (2013) 41 (D1): D94–D100. doi: 10.1093/nar/gks955

---

## Table of contents

1. [Quick start](#quick-start)
2. [Build](#build)
3. [Single-record run](#single-record-run)
4. [Multi-record FASTA: serial mode](#multi-record-fasta-serial-mode)
5. [Multi-record FASTA: parallel sharding](#multi-record-fasta-parallel-sharding)
6. [Output files](#output-files)
7. [Command-line options](#command-line-options)
8. [Performance notes](#performance-notes)
9. [Troubleshooting](#troubleshooting)

---

## Quick start

```sh
# 1. Build (uses gcc with -O2 -mssse3 by default; produces ./gfa)
make

# 2. Get the bundled test FASTA and golden outputs
tar -xf test_files.tar

# 3. Run on the example
./gfa -skipWGET -seq gfa_test.fasta -out gfa_test

# 4. The 14 output files (7 GFF + 7 TSV) should match test_files.tar exactly
for f in gfa_test_{IR,MR,DR,GQ,Z,STR,APR}.{gff,tsv}; do
    diff -q "$f" "$f"  # against your own — adjust path if you've separated golden
done
```

If `make` succeeds and the example run produces 14 files identical to the
ones inside `test_files.tar`, the installation is good.

---

## Build

`gfa` is a single C binary that depends only on libc and (at runtime) the
host CPU's SIMD ISA.

### Default build (SSSE3, runs on every x86_64 CPU since 2006)

```sh
make
```

Produces `./gfa`. The default `CFLAGS` in the Makefile are `-O2 -mssse3`.
SSSE3 has been standard on x86_64 since the Core 2 generation (2006), so
the resulting binary runs on any modern x86_64 machine.

### Optional: enable AVX2 path

```sh
make clean
make CFLAGS="-O2 -march=native"
```

This enables the AVX2 (32-byte) code paths in the SIMD inner loops. On
some Xeon-class CPUs the AVX clock-licensing penalty makes AVX2 slower
than SSSE3 for this workload, so we don't enable it by default. On
desktop / laptop CPUs without that penalty, `-march=native` is usually
~10–20 % faster.

### Optional: portable build (no SIMD, scalar fallbacks only)

```sh
make clean
make CFLAGS="-O2"
```

This disables the SIMD paths in `simd_match.h` and the
`cdna`/`rcdna` shuffle-table. The binary will be ~20× slower on the
inner-loop motif searches (IR/MR/DR), but it will compile and run on
non-x86 hosts.

### Cleaning up

```sh
make clean   # removes *.o; rerun `make` to rebuild
```

---

## Single-record run

The simplest invocation: one FASTA record in, seven `*.gff` and seven
`*.tsv` files out (one pair per motif type).

```sh
./gfa -skipWGET -seq INPUT.fasta -out OUTPUT_PREFIX
```

Required switches:
- `-seq <file>` — input FASTA file.
- `-out <prefix>` — output filename prefix. The seven motif suffixes
  (`_IR`, `_MR`, `_DR`, `_GQ`, `_Z`, `_STR`, `_APR`) and the file
  extensions (`.gff`, `.tsv`) are appended automatically.

Useful additional flags:
- `-skipWGET` — **almost always include this**. Without it, `gfa`
  attempts to `wget` an NCI URL at the end of the run to signal a web
  submission has completed. For local use this is just a wasted call.
- `-chrom <name>` — sequence name to use in column 1 of GFF/TSV. If
  omitted, the first whitespace-delimited token of the FASTA header is
  used.

Example:
```sh
./gfa -skipWGET -seq chr22.fa -out chr22 -chrom chr22
```

---

## Multi-record FASTA: serial mode

If the input FASTA contains multiple records, the default behaviour
processes all of them in sequence and **concatenates** their results into
each motif's GFF / TSV file:

```sh
./gfa -skipWGET -seq genome.fa -out genome
# -> genome_IR.gff contains rows from every record, in input order
# -> genome_IR.tsv likewise, with a single header line at the top
# -> ditto for MR / DR / GQ / Z / STR / APR
```

The per-record sequence name in column 1 comes from each FASTA header.
If `-chrom NAME` is given for a multi-record input, `gfa` appends
`_<record_index>` to disambiguate (e.g. `NAME_1`, `NAME_2`, …).

---

## Multi-record FASTA: parallel sharding

For multi-record inputs (whole genomes, multi-chromosome assemblies, etc.)
you can fan out one `gfa` worker per FASTA record, with each worker
processing a single record independently. This gives near-linear speedup
in the number of records up to the number of available cores.

Two equivalent ways to do this.

### Using the bundled wrapper (`gfa_parallel.sh`)

```sh
# -j N runs N workers in parallel; stitches the per-record outputs back
# into the requested -out prefix when done.
./gfa_parallel.sh -j 8 -seq genome.fa -out genome
```

What this does:
1. Counts records in `genome.fa`.
2. Launches up to `N` (`-j 8`) workers under `xargs -P`, each invoked as
   `./gfa -skipWGET -seq genome.fa -out genome_recK -record K`.
3. Concatenates `genome_recK_<TYPE>.{gff,tsv}` into `genome_<TYPE>.{gff,tsv}`
   for each motif TYPE, dropping duplicate `##gff-version` and TSV column
   headers so the stitched files are valid.
4. Removes the per-shard files when stitching succeeds.

Pass additional `gfa` options after `-out` and they are forwarded to every
worker. **Do not** pass `-record`, `-chrom`, or `-skipWGET` to the wrapper;
the wrapper manages those itself.

Forwarded options example:
```sh
./gfa_parallel.sh -j 16 -seq genome.fa -out genome -minDRrep 8 -maxDRrep 200
```

The wrapper looks for `./gfa` in the current directory. To use a different
binary, set the `GFA` environment variable:
```sh
GFA=/path/to/gfa ./gfa_parallel.sh -j 8 -seq genome.fa -out genome
```

### Rolling your own orchestrator (GNU parallel, Slurm, make, …)

`-record N` makes `gfa` process **only** the Nth FASTA record (1-based),
so any orchestrator that can iterate `1..N` can drive a fan-out. Each
worker must use a **unique `-out` prefix** so they don't collide on the
same output file.

GNU parallel:
```sh
NREC=$(grep -c '^>' genome.fa)
seq 1 "$NREC" | parallel -j 8 \
    ./gfa -skipWGET -seq genome.fa -out genome_rec{} -record {}
```

Slurm array job:
```sh
#SBATCH --array=1-93
#SBATCH --cpus-per-task=1
./gfa -skipWGET -seq genome.fa -out genome_rec${SLURM_ARRAY_TASK_ID} \
      -record ${SLURM_ARRAY_TASK_ID}
```

`make -j`:
```makefile
NREC := $(shell grep -c '^>' genome.fa)
SHARDS := $(addprefix genome_rec,$(shell seq 1 $(NREC)))

all: $(SHARDS)

genome_rec%_IR.gff: genome.fa
	./gfa -skipWGET -seq $< -out genome_rec$* -record $*
```

After all workers complete, concatenate per-type outputs in the order you
prefer (the wrapper script's stitching loop is a good template).

---

## Output files

For each motif type the run produces two files:

| Suffix | Content |
|---|---|
| `_IR.gff` / `_IR.tsv` | Inverted Repeats (cruciform-forming) |
| `_MR.gff` / `_MR.tsv` | Mirror Repeats (triplex-forming) |
| `_DR.gff` / `_DR.tsv` | Direct Repeats (slipped DNA) |
| `_GQ.gff` / `_GQ.tsv` | G-Quadruplexes |
| `_Z.gff` / `_Z.tsv` | Z-DNA |
| `_STR.gff` / `_STR.tsv` | Short Tandem Repeats |
| `_APR.gff` / `_APR.tsv` | A-Phased Repeats (bent DNA) |

GFF files follow [GFF2](http://gmod.org/wiki/GFF2) (9 tab-separated columns
plus a `##gff-version 3` header line). TSV files have a single header row
followed by one motif per line. Coordinates are 1-based, inclusive.

You can suppress any motif's search with the matching `-skip<TYPE>` flag —
the output files for that type are simply not produced.

---

## Command-line options

### Required

| Switch | Description |
|---|---|
| `-seq <file>` | Input DNA FASTA file. |
| `-out <prefix>` | Output filename prefix. Motif abbreviations and file extensions are appended automatically. |

### Sharding

| Switch | Description |
|---|---|
| `-record <N>` | Process only the Nth FASTA record (1-based). Use a unique `-out` prefix per shard so concurrent invocations do not collide. Intended for node-scope parallelism via `gfa_parallel.sh` or an external orchestrator. |

### Numeric thresholds (each takes a value; defaults shown)

```
-minGQrep <3>          min consecutive G's per G run                (no max)
-maxGQspacer <7>       max distance between G runs                  (min 1)
-minMRrep <10>         min length of half a mirror repeat           (no max)
-maxMRspacer <100>     max distance between mirror repeat halves    (min 0)
-minIRrep <6>          min length of half an inverted repeat        (no max)
-maxIRspacer <100>     max distance between inverted repeat halves  (min 0)
-shortIRcut <9>        max length to classify an IR as "short"
-shortIRspacer <4>     max distance for short inverted repeats      (min 0)
-minDRrep <10>         min length of half a direct repeat
-maxDRrep <300>        max length of half a direct repeat
-maxDRspacer <10>      max distance between direct repeat halves    (min 0)
-minATracts <3>        min consecutive A tracts to form an APR
-minATractSep <10>     min separation between A tract centres
-maxATractSep <11>     max separation between A tract centres
-maxAPRlen <9>         max consecutive A's in an A tract
-minAPRlen <3>         min consecutive A's in an A tract
-minZlen <10>          min length of Z-DNA pur/pyr run              (no max)
-minSTR <1>            min length of STR repeating element
-maxSTR <9>            max length of STR repeating element
-minSTRbp <8>          min overall length for STR qualification
-minCruciformRep <6>   min IR repeat length to qualify as cruciform
-maxCruciformSpacer <4> max IR spacer to qualify as cruciform
-minTriplexYRpercent <10> min pur/pyr % for MR to qualify as triplex
-maxTriplexSpacer <8>  max MR spacer to qualify as triplex
-maxSlippedSpacer <0>  max DR spacer to qualify as slipped
```

### Boolean / value-less switches

```
-chrom <name>     Sequence name for column 1 of GFF/TSV
                  (default: first whitespace-delimited token of fasta header)
-skipAPR          Don't search for A-phased repeats (bent DNA)
-skipSTR          Don't search for short tandem repeats
-skipDR           Don't search for direct repeats (slipped DNA)
-skipMR           Don't search for mirror repeats (triplex DNA)
-skipIR           Don't search for inverted repeats (cruciform DNA)
-skipGQ           Don't search for G-quadruplex motifs
-skipZ            Don't search for Z-DNA motifs
-skipSlipped      Don't compute the slipped subset of DRs
-skipCruciform    Don't compute the cruciform subset of IRs
-skipTriplex      Don't compute the triplex subset of MRs
-skipWGET         Don't make the wget call to the NCI PHP endpoint
                  (recommended for local use)
-doCHMOD          Run `chmod 664` on the output files after writing
```

Print the full list at runtime with `./gfa` (no arguments).

---

## Performance notes

The current binary is heavily optimised for x86_64 with SSSE3 SIMD:

- DNA, REP and island buffers are heap-allocated and sized to the actual
  longest record (no static 300 Mbp / 1.5 GB BSS cap).
- The FASTA reader uses buffered `fread` + a 256-byte lookup table
  (no `getc` per byte).
- The IR / MR / DR motif finders precompute a per-record LCP table from
  16-D-wide SIMD broadcast-match chunks and short-circuit when no LCP is
  large enough to satisfy the motif minimum. On uniform random DNA the
  short-circuit fires for ~98–99.97 % of positions.
- STR uses one SIMD `forward_match` call per (i, rpsz) pair instead of a
  `strncmp` loop, deriving the rep count as `1 + LCP/rpsz`.
- `cdna` / `rcdna` use a `_mm_shuffle_epi8` LUT keyed on the low nibble
  of each ASCII base, complementing 16 bytes per instruction.

Measured on a Xeon @ 2.8 GHz against the pre-optimisation baseline:

| Input | Before | After | Speedup |
|---|---|---|---|
| 1 MB random FASTA | 11.15 s | 0.51 s | **22 ×** |
| 10 MB random FASTA | ~111 s (extrapolated) | 5.1 s | **22 ×** |
| 1 MB synthetic FASTA with planted inverted repeats | 3.13 s | 0.12 s | **26 ×** |

Memory: peak RSS dropped from ~900 MB to ~2 MB on the bundled test input,
and the hard 300 Mbp input cap is gone.

For multi-record inputs, `gfa_parallel.sh` multiplies these wins by the
number of cores you give it.

---

## Troubleshooting

**`make` fails with "no SSSE3" or similar.** You're on a CPU older than
Core 2 (2006), or in a build environment that overrides the compiler flags.
Try the portable build: `make clean && make CFLAGS="-O2"`.

**`Illegal instruction` (SIGILL) at runtime.** The binary was built with an
ISA the host doesn't support — most commonly, built with
`-march=native` on a machine with AVX2 and then moved to a host without
AVX2. Rebuild on (or for) the target host, or use the default
`-mssse3` build which runs everywhere.

**`FATAL ERROR: Sequence file NOT in FASTA format`.** The first byte of
the file isn't `>`. Strip any BOM / leading whitespace and ensure the file
really is a FASTA.

**`-record N exceeds number of FASTA records (M)`.** You passed `-record`
with an index past the end of the input. Note that `-record` is 1-based,
not 0-based.

**Sharded outputs from `gfa_parallel.sh` look truncated.** A worker probably
errored. The wrapper writes per-shard stderr to `${OUT}_rec<N>.stderr` and
removes those files only on a clean run; if they remain, inspect them for
the failure.

---

## Authors and history

- Original author: Duncan E. Donohue, Ph.D.
- Earlier work: Jack R. Collins, Ph.D.
- Performance / parallelism work: see `git log` on this branch.
