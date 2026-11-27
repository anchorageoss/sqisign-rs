#!/bin/sh
# Run every check of paper-shortiso and write one log per item to
# paper-shortiso/results/. Exit status 0 only if every test passed.
#
#   sh paper-shortiso/run.sh
#
# Records (class polynomials, endpoint orders) are compared with the files
# checked in; SHORTISO_WRITE=1 in the environment regenerates them.
set -eu
cd "$(dirname "$0")/.."
out=paper-shortiso/results
mkdir -p "$out"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
status=0

cargo test --release -p sqisign-rs --no-run \
    --test shortiso_sqisign --test shortiso_generic \
    --test shortiso_equations --test shortiso_table \
    --test shortiso_minimum >/dev/null 2>&1

# run one test binary (further arguments go to the test harness); keep its
# output from "running N tests" on
run() {
    name=$1
    shift
    if ! cargo test --release -p sqisign-rs --test "$name" -- --nocapture "$@" \
        >"$tmp/$name.raw" 2>&1; then
        status=1
    fi
    sed -n '/^running [0-9]* tests*$/,$p' "$tmp/$name.raw" | sed 's/|[0-9.]*s$//' >"$tmp/$name.log"
}

header() {
    echo "# paper-shortiso, $1"
    echo "# $(git rev-parse --short HEAD) $(date -u +%Y-%m-%dT%H:%MZ) $(rustc --version)"
    echo "# produced by paper-shortiso/run.sh from $2"
}

run shortiso_sqisign
run shortiso_generic
run shortiso_equations
run shortiso_table
run shortiso_minimum --test-threads=1

{
    header "item 1: the three SQIsign primes on the crate's field, curve, basis and stored action" crates/sqisign-rs/tests/shortiso_sqisign.rs
    grep -E '^ROW\|(item1|item3-cross)\||^test |panicked|^assertion|^ +(left|right):' "$tmp/shortiso_sqisign.log" || true
} >"$out/item1.log"
{
    header "item 2: the quaternion side at the three SQIsign primes" crates/sqisign-rs/tests/shortiso_sqisign.rs
    grep -E '^ROW\|item2\||^test |panicked|^assertion|^ +(left|right):' "$tmp/shortiso_sqisign.log" || true
} >"$out/item2.log"
{
    header "item 3: the other primes over the generic field" crates/sqisign-rs/tests/shortiso_generic.rs
    cat "$tmp/shortiso_generic.log"
} >"$out/item3.log"
{
    header "item 4: the norm equations, integers only" crates/sqisign-rs/tests/shortiso_equations.rs
    cat "$tmp/shortiso_equations.log"
} >"$out/item4.log"
{
    header "item 5: what stays outside Rust" "a search of the repository (no Rust in this item)"
    echo "Sage scripts of the paper searched for under paper-shortiso/ and the repository:"
    found=$(git ls-files | grep -i -E 'shortiso.*\.(sage|py)$|shortest.*\.(sage|py)$' || true)
    if [ -n "$found" ]; then
        echo "$found"
    else
        echo "none found, none attached to the request. Nothing was run in Sage or Python."
    fi
    echo "Sage only, not run here: check N1 (least degrees from classical modular polynomials,"
    echo "independent of the endomorphism ring), and the kernel check at small primes whose"
    echo "torsion needs an extension of F_p of degree above 2."
    echo "See the coverage table in paper-shortiso/REPORT.md."
} >"$out/item5.log"
{
    header "item 6: the (D, m, h(-D)) table and the claims about other schemes" crates/sqisign-rs/tests/shortiso_table.rs
    cat "$tmp/shortiso_table.log"
    echo
    echo "Source references for the claims about SQIsign and MIKE: paper-shortiso/SOURCES.md"
} >"$out/item6.log"

{
    header "item 7: the minimality of Corollary A by exact enumeration in the ideal classes" crates/sqisign-rs/tests/shortiso_minimum.rs
    cat "$tmp/shortiso_minimum.log"
} >"$out/item7.log"

# the table of primes
{
    printf 'group\tprime\tbits\tD\tm\tsteps\tendpoint\tchain\tresult\n'
    grep -h '^ROW|item1|' "$tmp/shortiso_sqisign.log" |
        awk -F'|' '{ sub(/^D=/, "", $5); sub(/^m=/, "", $6); n = $6; sub(/^2\^/, "", n);
            printf "sqisign, crate machinery\t%s\t-\t%s\t%s\t%s\t%s\t%s\t%s\n", $4, $5, $6, n, $8, $9, $NF }'
    grep -h '^ROW|' "$tmp/shortiso_generic.log" |
        awk -F'|' '{ sub(/ bits$/, "", $4); sub(/^D=/, "", $5); sub(/^m=/, "", $6); sub(/^steps=/, "", $7);
            printf "%s, generic field\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $5, $6, $7, $8, $9, $10 }' |
        sort -t "$(printf '\t')" -k1,1 -k3,3n
} >"$out/primes.tsv"

for name in shortiso_sqisign shortiso_generic shortiso_equations shortiso_table shortiso_minimum; do
    if ! grep -q '^test result: ok' "$tmp/$name.log"; then
        status=1
        echo "FAILED: $name (see $out)" >&2
        tail -n 40 "$tmp/$name.raw" >&2
    fi
done
[ "$status" -eq 0 ] && echo "paper-shortiso: all checks passed; logs in $out"
exit "$status"
