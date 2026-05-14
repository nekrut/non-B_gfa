#!/bin/sh
# gfa_parallel.sh -- run gfa on each FASTA record in parallel, then stitch outputs.
#
# Usage:  gfa_parallel.sh -j <jobs> -seq <fasta> -out <prefix> [other gfa args...]
#
# Fans out one gfa worker per FASTA record using `-record N`. Per-shard outputs
# are written to <prefix>_rec<N>_<TYPE>.{gff,tsv} and then concatenated into
# <prefix>_<TYPE>.{gff,tsv}. Shard files are removed on success.
#
# Pass any additional gfa switches after -out; they are forwarded to every
# worker. Do NOT pass -record, -chrom (per-record chrom is derived from the
# fasta title), or -skipWGET (it is forced on to avoid N concurrent wget calls).

set -e

JOBS=4
SEQ=
OUT=
GFA="${GFA:-./gfa}"
PASSTHROUGH=

while [ $# -gt 0 ]; do
    case "$1" in
        -j)   JOBS="$2"; shift 2 ;;
        -seq) SEQ="$2"; shift 2 ;;
        -out) OUT="$2"; shift 2 ;;
        -record|-chrom|-skipWGET)
            echo "$0: $1 is managed by the wrapper and cannot be passed through" >&2
            exit 2 ;;
        *)    PASSTHROUGH="$PASSTHROUGH $1"; shift ;;
    esac
done

if [ -z "$SEQ" ] || [ -z "$OUT" ]; then
    echo "Usage: $0 -j <jobs> -seq <fasta> -out <prefix> [other gfa args]" >&2
    exit 2
fi
if [ ! -x "$GFA" ]; then
    echo "$0: $GFA not executable (set \$GFA to override)" >&2
    exit 2
fi
if [ ! -f "$SEQ" ]; then
    echo "$0: $SEQ not found" >&2
    exit 2
fi

NREC=$(grep -c '^>' "$SEQ")
if [ "$NREC" -lt 1 ]; then
    echo "$0: no FASTA records found in $SEQ" >&2
    exit 2
fi

echo "$0: $NREC records, $JOBS parallel workers" >&2

# Fan out: one worker per record. xargs -P limits concurrency.
seq 1 "$NREC" | xargs -n1 -P"$JOBS" -I{} sh -c \
    "$GFA -skipWGET -seq \"$SEQ\" -out \"${OUT}_rec{}\" -record {} $PASSTHROUGH \
        >\"${OUT}_rec{}.stdout\" 2>\"${OUT}_rec{}.stderr\""

# Stitch: concatenate per-type files. Order of records is preserved.
for TYPE in IR MR DR GQ Z STR APR; do
    for EXT in gff tsv; do
        OUTFILE="${OUT}_${TYPE}.${EXT}"
        : > "$OUTFILE"
        first=1
        for N in $(seq 1 "$NREC"); do
            SHARD="${OUT}_rec${N}_${TYPE}.${EXT}"
            [ -f "$SHARD" ] || continue
            if [ $first -eq 0 ]; then
                if [ "$EXT" = "gff" ]; then
                    # Skip the ##gff-version header on subsequent shards
                    grep -v '^##gff-version' "$SHARD" >> "$OUTFILE"
                else
                    # Skip the column header line on subsequent shards
                    tail -n +2 "$SHARD" >> "$OUTFILE"
                fi
            else
                cat "$SHARD" >> "$OUTFILE"
            fi
            first=0
            rm -f "$SHARD"
        done
        # Remove the empty stitched file if no shard produced output for this type
        [ -s "$OUTFILE" ] || rm -f "$OUTFILE"
    done
done

# Clean up worker logs
rm -f "${OUT}"_rec*.stdout "${OUT}"_rec*.stderr

echo "$0: done. Stitched outputs: ${OUT}_<TYPE>.{gff,tsv}" >&2
