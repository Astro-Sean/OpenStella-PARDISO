#!/usr/bin/env bash
# Standalone runner for a shell+wind model produced by
# make_shell_wind_model.py.  Chains the three stages STELLA needs:
#
#   eve/      eve2.exe      .hyd + .abn -> .mod + .xni (STELLA binary model)
#   vladsf/   xronfict.exe  .mod + .xni -> .1-.6 + .ab (opacity tables)
#   strad/    xstella.exe   .dat + .mod + .xni + .1 -> .swd .tt .res ...
#
# Usage: bash run_model.sh <workdir> [threads]
#   workdir  directory written by make_shell_wind_model.py
#   threads  OMP/MKL thread count (default 2)

set -euo pipefail

WORKDIR="$(cd "${1:?usage: run_model.sh <workdir> [threads]}" && pwd)"
NTHREADS="${2:-2}"

# STELLA_ROOT = the repo root (two levels up from this script).
STELLA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MODEL="$(basename "$(ls "$WORKDIR"/modmake/*.hyd 2>/dev/null | head -1)" .hyd)"
if [[ -z "$MODEL" ]]; then
    echo "ERROR: no .hyd found in $WORKDIR/modmake" >&2
    exit 1
fi

export HOMEStella="$STELLA_ROOT"
export OMP_NUM_THREADS="$NTHREADS"
export MKL_NUM_THREADS="$NTHREADS"
export MKL_THREADING_LAYER=intel
if [[ -n "${MKLROOT:-}" && -d "$MKLROOT/lib/intel64" ]]; then
    export LD_LIBRARY_PATH="$MKLROOT/lib/intel64:${LD_LIBRARY_PATH:-}"
fi

EVE_EXE="$STELLA_ROOT/run/eve/eve2.exe"
RONF_EXE="$STELLA_ROOT/run/vladsf/xronfict.exe"
XSTELLA="$STELLA_ROOT/bin/xstella.exe"
for exe in "$EVE_EXE" "$RONF_EXE" "$XSTELLA"; do
    [[ -x "$exe" ]] || { echo "ERROR: missing $exe - run build_parallel.sh first" >&2; exit 1; }
done

echo "======================================"
echo "Model: $MODEL   workdir: $WORKDIR"
echo "Threads: $NTHREADS   STELLA_ROOT: $STELLA_ROOT"
echo "======================================"

# --- Stage 1: eve2 (.hyd/.abn -> .mod/.xni) -------------------------------
echo ""
echo "--- Stage 1: eve2 ---"
cd "$WORKDIR/eve"
"$EVE_EXE" 2>&1 | tee eve2.log
[[ -s "$MODEL.mod" && -s "$MODEL.xni" ]] || {
    echo "ERROR: eve2 did not produce $MODEL.mod/.xni" >&2; exit 1; }

# --- Stage 2: xronfict (.mod/.xni -> opacity .1-.6 + .ab) ------------------
echo ""
echo "--- Stage 2: xronfict (opacity tables, ~15-20 min) ---"
cd "$WORKDIR/vladsf"
ln -sfn "$STELLA_ROOT/vladsf/lineatom.dat" lineatom.dat
ln -sfn "$STELLA_ROOT/run/vladsf/yakovlev" yakovlev
# A pre-compiled dump skips the line-list parse; only valid for the
# default lineatom.dat.
DUMP="$STELLA_ROOT/models/X_STELA/vladsf/linedata.dump"
[[ -f "$DUMP" ]] && ln -sf "$DUMP" linedata.dump
"$RONF_EXE" 2>&1 | tee ronfict.log
for nf in 1 2 3 4 5 6; do
    [[ -s "$MODEL.$nf" ]] || {
        echo "ERROR: missing opacity file $MODEL.$nf" >&2; exit 1; }
done
[[ -s "$MODEL.ab" ]] || { echo "ERROR: missing $MODEL.ab" >&2; exit 1; }

# --- Stage 3: xstella ------------------------------------------------------
echo ""
echo "--- Stage 3: xstella ---"
cd "$WORKDIR/strad"
cp "$WORKDIR/modmake/$MODEL.dat" .
ln -sf "../eve/$MODEL.mod" "$MODEL.mod"
ln -sf "../eve/$MODEL.xni" "$MODEL.xni"
ln -sf "../vladsf/$MODEL.1" "$MODEL.1"
# STELLA cannot overwrite outputs; clear leftovers unless restarting.
rm -f "$MODEL".{res,bm,balance,swd,swd_rad,tt,flx,lbol,ph,prf,crv}
"$XSTELLA" < /dev/null 2>&1 | tee stella_run.log

# End-of-run outputs are moved to $HOMEStella/res/; pull them back here.
for ext in swd swd_rad tt lbol ph flx res; do
    cp "$STELLA_ROOT/res/$MODEL.$ext" "$WORKDIR/strad/" 2>/dev/null || true
done

echo ""
echo "======================================"
echo "Done.  Outputs in $WORKDIR/strad/ and $STELLA_ROOT/res/"
echo "  $MODEL.swd  structure + composition per epoch"
echo "  $MODEL.tt   15-band light curve"
echo "======================================"
