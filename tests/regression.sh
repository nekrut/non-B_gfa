#!/usr/bin/env bash
# tests/regression.sh -- non-B-gfa regression test suite.
#
# Run from the project root (so that ./gfa is built and test_files.tar exists)
# OR from anywhere with the repo as the first argument.
#
# Exits 0 on success, non-zero on any failure.

set -e

# Allow running from anywhere: cd into the project root.
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

GFA="${GFA:-./gfa}"
if [ ! -x "$GFA" ]; then
    echo "$0: $GFA not found or not executable. Run \`make\` first." >&2
    exit 2
fi

if [ ! -f test_files.tar ]; then
    echo "$0: test_files.tar not found in $(pwd)" >&2
    exit 2
fi

# Use a temp work area so the project root stays clean.
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

PASS=0
FAIL=0

fail() {
    echo "  FAIL: $1" >&2
    FAIL=$((FAIL + 1))
}

ok() {
    echo "  ok:   $1"
    PASS=$((PASS + 1))
}

# ----------------------------------------------------------------------------
# Test 1: bundled test_files.tar golden
#
# Run gfa with defaults on gfa_test.fasta; outputs must be byte-for-byte
# identical to the seven .gff + seven .tsv files inside test_files.tar.
# ----------------------------------------------------------------------------

echo "Test 1: bundled gfa_test.fasta golden"

mkdir -p "$WORK/t1/golden" "$WORK/t1/got"
tar -xf test_files.tar -C "$WORK/t1/golden"
cp "$WORK/t1/golden/gfa_test.fasta" "$WORK/t1/got/"

( cd "$WORK/t1/got" && "$REPO_ROOT/$GFA" -skipWGET -seq gfa_test.fasta -out gfa_test ) \
    >"$WORK/t1.stdout" 2>"$WORK/t1.stderr"

for f in gfa_test_IR.gff  gfa_test_IR.tsv  \
         gfa_test_MR.gff  gfa_test_MR.tsv  \
         gfa_test_DR.gff  gfa_test_DR.tsv  \
         gfa_test_GQ.gff  gfa_test_GQ.tsv  \
         gfa_test_Z.gff   gfa_test_Z.tsv   \
         gfa_test_STR.gff gfa_test_STR.tsv \
         gfa_test_APR.gff gfa_test_APR.tsv; do
    if diff -q "$WORK/t1/got/$f" "$WORK/t1/golden/$f" >/dev/null 2>&1; then
        ok "golden $f"
    else
        fail "golden $f differs from test_files.tar"
    fi
done

# ----------------------------------------------------------------------------
# Test 2: -record sharded vs serial multi-record run
#
# Build a synthetic 3-record fasta by concatenating the test fasta three
# times under different headers, then verify gfa_parallel.sh produces the
# same per-motif lines (canonicalised by sort) as a single serial gfa run.
# ----------------------------------------------------------------------------

echo "Test 2: serial vs sharded multi-record run"

mkdir -p "$WORK/t2/serial" "$WORK/t2/parallel"
awk 'BEGIN{n=0} /^>/{n++; print ">rec"n; next} {print}' \
    "$WORK/t1/golden/gfa_test.fasta" >  "$WORK/t2/rec123.fa"
awk 'BEGIN{n=1} /^>/{n++; print ">rec"n; next} {print}' \
    "$WORK/t1/golden/gfa_test.fasta" >> "$WORK/t2/rec123.fa"
awk 'BEGIN{n=2} /^>/{n++; print ">rec"n; next} {print}' \
    "$WORK/t1/golden/gfa_test.fasta" >> "$WORK/t2/rec123.fa"

( cd "$WORK/t2" && "$REPO_ROOT/$GFA" -skipWGET -seq rec123.fa -out serial/all ) \
    >"$WORK/t2-serial.stdout" 2>"$WORK/t2-serial.stderr"

( cd "$WORK/t2" && GFA="$REPO_ROOT/$GFA" "$REPO_ROOT/gfa_parallel.sh" \
        -j 3 -seq rec123.fa -out parallel/all ) \
    >"$WORK/t2-parallel.stdout" 2>"$WORK/t2-parallel.stderr"

for f in "$WORK/t2/serial"/all_*.tsv; do
    base=$(basename "$f")
    if diff -q \
            <(sort "$f") \
            <(sort "$WORK/t2/parallel/$base") >/dev/null 2>&1; then
        ok "sharded $base"
    else
        fail "sharded $base differs from serial"
    fi
done

for f in "$WORK/t2/serial"/all_*.gff; do
    base=$(basename "$f")
    if diff -q \
            <(grep -v '^##gff-version' "$f"                 | sort) \
            <(grep -v '^##gff-version' "$WORK/t2/parallel/$base" | sort) \
            >/dev/null 2>&1; then
        ok "sharded $base"
    else
        fail "sharded $base differs from serial"
    fi
done

# ----------------------------------------------------------------------------
# Test 3: -record argument validation
#
# `-record -5` and `-record 0` must both be rejected with a FATAL ERROR.
# ----------------------------------------------------------------------------

echo "Test 3: -record validation"

for bad in -5 0; do
    if "$GFA" -skipWGET -seq "$WORK/t1/golden/gfa_test.fasta" \
            -out "$WORK/t3" -record "$bad" >/dev/null 2>"$WORK/t3.stderr"; then
        fail "-record $bad was accepted (should be rejected)"
    elif grep -q "FATAL ERROR" "$WORK/t3.stderr"; then
        ok "-record $bad rejected"
    else
        fail "-record $bad rejected but without FATAL ERROR (got: $(head -1 "$WORK/t3.stderr"))"
    fi
done

# ----------------------------------------------------------------------------
# Test 4: IUPAC ambiguity-code handling
#
# Build an input containing q/s/d/w (the IUPAC codes whose ASCII low nibble
# collides with a/c/t/g) and verify that the result is identical to a build
# with the SIMD complement path disabled. Skipped if we can't build a
# secondary scalar binary -- e.g. on CI without a separate gcc.
# ----------------------------------------------------------------------------

echo "Test 4: IUPAC ambiguity-code handling"

cat > "$WORK/t4-iupac.fa" <<'EOF'
>iupac_mix
acgtnswdqacgtnswdqacgtnswdqacgtnswdqacgtnswdqacgtnswdqacgtnswdqac
acgtnacgtnacgtnacgtnacgtnacgtnacgtnacgtnacgtnacgtnacgtnacgtnacgtn
EOF

# Build a scalar (no -mssse3) gfa next to the one under test so we can diff
# their outputs.
SCALAR="$WORK/gfa_scalar"
(
    cd "$WORK"
    cp "$REPO_ROOT"/*.c "$REPO_ROOT"/*.h "$REPO_ROOT/Makefile" .
    make CFLAGS="-O2" clean >/dev/null 2>&1 || true
    make CFLAGS="-O2" >"$WORK/t4-build.log" 2>&1
) || true

if [ -x "$WORK/gfa" ]; then
    mv "$WORK/gfa" "$SCALAR"
    mkdir -p "$WORK/t4/scalar" "$WORK/t4/simd"
    "$SCALAR" -skipWGET -seq "$WORK/t4-iupac.fa" -out "$WORK/t4/scalar/r" \
        >/dev/null 2>&1
    "$GFA"    -skipWGET -seq "$WORK/t4-iupac.fa" -out "$WORK/t4/simd/r"   \
        >/dev/null 2>&1

    for f in "$WORK/t4/scalar"/r_*.gff "$WORK/t4/scalar"/r_*.tsv; do
        [ -f "$f" ] || continue
        base=$(basename "$f")
        if diff -q "$f" "$WORK/t4/simd/$base" >/dev/null 2>&1; then
            ok "iupac $base (SIMD == scalar)"
        else
            fail "iupac $base SIMD output differs from scalar"
        fi
    done
else
    echo "  skip: could not build scalar gfa for IUPAC comparison" >&2
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------

echo
echo "=========================================="
echo "regression.sh: $PASS passed, $FAIL failed"
echo "=========================================="

if [ "$FAIL" -ne 0 ]; then
    exit 1
fi
exit 0
