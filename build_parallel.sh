#!/usr/bin/env bash
# Build script for OpenStella_parallel (ifort + MKL PARDISO + OpenMP).
# Uses pre-generated .f files (no trefor preprocessing needed for the solver path).
#
# Usage:  bash build_parallel.sh <NZON>
#   NZON  = radial zone count of your model (sets Mzon). REQUIRED for correctness;
#           defaults to 200, which under-dimensions larger models -> pass it.
#
# Requirements: Intel Fortran (ifort or ifx) and Intel MKL (MKLROOT set).
#
# The script detects the host and sets up the toolchain automatically:
#   viper* / raven*  -> load the site Intel module (ifort + MKL)
#   anything else    -> source the standard oneAPI setvars.sh if needed
# Override detection with STELLA_SITE=viper|raven|generic.

set -uo pipefail

# This script lives at the root of the OpenStella_parallel tree.
STELLA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NZON="${1:-200}"

echo "======================================"
echo "Building OpenStella_parallel (Intel Fortran + MKL + OpenMP, NZON=$NZON)"
echo "STELLA_ROOT: $STELLA_ROOT"
echo "======================================"

if [[ $# -lt 1 ]]; then
    echo "WARNING: NZON not given; defaulting to Mzon=200."
    echo "         If your model uses more than 200 radial zones, re-run with:"
    echo "         bash build_parallel.sh <NZON>"
fi

# --- Step 0: compiler + MKL environment ------------------------------------
# MPCDF login nodes are named viperNN / ravenNN. On those, the Intel module
# provides ifort and MKLROOT in one step. Elsewhere we rely on the oneAPI
# setvars.sh that a normal install drops under /opt/intel or $HOME.
HOST_SHORT="$(hostname -s 2>/dev/null || hostname)"
SITE="${STELLA_SITE:-}"
if [[ -z "$SITE" ]]; then
    case "$HOST_SHORT" in
        viper*)      SITE=viper ;;
        raven*|rav*) SITE=raven ;;
        *)           SITE=generic ;;
    esac
fi
echo "Host: $HOST_SHORT (site: $SITE)"

load_intel_module() {
    command -v module >/dev/null 2>&1 || return 1
    local mod
    for mod in "$@"; do
        if module load "$mod" 2>/dev/null; then
            echo "Loaded module: $mod"
            return 0
        fi
    done
    return 1
}

case "$SITE" in
    viper|raven)
        # Same MPCDF module tree on both machines; try newest-compatible
        # first, then the generic names.
        load_intel_module intel/2024.0 intel-oneapi intel || true
        ;;
    *)
        if [[ -z "${MKLROOT:-}" ]]; then
            for vars in \
                /opt/intel/oneapi/setvars.sh \
                "$HOME/intel/oneapi/setvars.sh" \
                "$HOME/.intel/oneapi/setvars.sh"; do
                if [[ -f "$vars" ]]; then
                    echo "Sourcing $vars"
                    # setvars.sh prints banners; keep the build log clean.
                    # shellcheck disable=SC1090
                    source "$vars" >/dev/null 2>&1 || true
                    break
                fi
            done
        fi
        ;;
esac

# --- sanity checks ---------------------------------------------------------
# ifx is the LLVM-based replacement shipped since oneAPI 2024; the same
# flags below work for both drivers.
if command -v ifort >/dev/null 2>&1; then
    FC=ifort
elif command -v ifx >/dev/null 2>&1; then
    FC=ifx
else
    echo "ERROR: no Intel Fortran compiler (ifort/ifx) in PATH." >&2
    echo "  MPCDF:    module load intel/2024.0" >&2
    echo "  other PC: source /opt/intel/oneapi/setvars.sh" >&2
    exit 1
fi
echo "FC: $FC"
if [[ -z "${MKLROOT:-}" ]]; then
    echo "ERROR: MKLROOT not set. Source the Intel MKL environment (setvars.sh)." >&2
    exit 1
fi
echo "MKLROOT: $MKLROOT"

# Keep the MKL runtime on LD_LIBRARY_PATH for the linked executable; a
# module purge elsewhere may not have restored it.
if [[ -d "$MKLROOT/lib/intel64" ]]; then
    case ":${LD_LIBRARY_PATH:-}:" in
        *":$MKLROOT/lib/intel64:"*) ;;
        *) export LD_LIBRARY_PATH="$MKLROOT/lib/intel64:${LD_LIBRARY_PATH:-}" ;;
    esac
fi

export SYSTYPE="ifort"
export HOMEStella="$STELLA_ROOT"
mkdir -p "$STELLA_ROOT/bin"
export PATH="$STELLA_ROOT/bin:$PATH"

# --- Step 1: zone.inc ------------------------------------------------------
echo ""
echo "--- Step 1: Set zone.inc (Mzon=$NZON) ---"
cd "$STELLA_ROOT/src"
if [[ ! -f zone.inc ]]; then
    cp zoneSample.inc zone.inc
fi
cp opacityHomo.inc opacity.inc
sed -i "s/Mzon=[0-9]*/Mzon=$NZON/g" zone.inc
echo "zone.inc Mzon line:"; grep "Mzon" zone.inc | head -2

# --- Step 2: compile -------------------------------------------------------
echo ""
echo "--- Step 2: Compile STELLA objects ---"
cd "$STELLA_ROOT/obj"

# Symlink include files so the compiler finds them.
for f in "$STELLA_ROOT"/src/*.inc; do
    bname=$(basename "$f")
    [[ ! -e "$bname" ]] && ln -sf "$f" "$bname"
done
for f in "$STELLA_ROOT"/vladsf/*.inc; do
    bname=$(basename "$f")
    [[ ! -e "$bname" ]] && ln -sf "$f" "$bname" 2>/dev/null || true
done

FFLAGS="-c -132 -save -zero -O3 -fp-model strict -qopenmp -D__INTEL"
FFLAGS_F90="-c -save -zero -O3 -fp-model strict -qopenmp -D__INTEL"
INCL_DIR="-I$STELLA_ROOT/src/ -I$STELLA_ROOT/vladsf -I${MKLROOT}/include -I$STELLA_ROOT/obj"
MKL_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl"

SOURCES="stradsep5tt.f begradsep.f cosetbgh.f hcdfjrad.f \
        hcdfnrad.f traneq.f eddi.f gdepos6.f nthnew.f \
        stiffbghY12m.f lbalsw.f stradio.f \
        vtimef90.f sahaandd.f ubv.f obsubvri.f \
        tt4strad.f begtt.f lbol.f \
        burnc.f volenpumnoint.f hapsepnc.f hcdhaph.f \
        oparon.f length.f words.f azdat.f pardiso_wrap.f90"

rm -f *.o ../run/strad/xstella6new.exe

find_source() {
    local fname="$1"
    for dir in "$STELLA_ROOT/strad" "$STELLA_ROOT/src" "$STELLA_ROOT/vladsf"; do
        if [[ -f "$dir/$fname" ]]; then
            echo "$dir/$fname"
            return
        fi
    done
    echo ""
}

COMPILE_ERRORS=0

# Compile the PARDISO module first so its .mod is available to the integrator.
echo "Compiling pardiso_wrap.f90 (module)..."
PARDISO_SRC="$STELLA_ROOT/strad/pardiso_wrap.f90"
if [[ -f "$PARDISO_SRC" ]]; then
    if $FC $FFLAGS_F90 "$PARDISO_SRC" $INCL_DIR -o pardiso_wrap.o 2>&1 | grep -v "remark #10448"; then
        :
    else
        echo "ERROR compiling pardiso_wrap.f90" >&2
        COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
    fi
else
    echo "ERROR: Cannot find pardiso_wrap.f90" >&2
    COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
fi

for src in $SOURCES; do
    if [[ "$src" == "pardiso_wrap.f90" ]]; then
        continue
    fi
    obj="${src%.f}.o"
    srcpath=$(find_source "$src")
    if [[ -z "$srcpath" ]]; then
        echo "ERROR: Cannot find $src" >&2
        COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
        continue
    fi
    echo "Compiling $src..."
    if $FC $FFLAGS "$srcpath" $INCL_DIR -o "$obj" 2>&1 | grep -v "remark #10448"; then
        :
    else
        echo "ERROR compiling $src" >&2
        COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
    fi
done

if [[ $COMPILE_ERRORS -gt 0 ]]; then
    echo "ERROR: $COMPILE_ERRORS compilation failures" >&2
    exit 1
fi
echo "All objects compiled."

# --- Step 3: link ----------------------------------------------------------
echo ""
echo "--- Step 3: Link xstella6new.exe ---"
OBJS=""
for src in $SOURCES; do
    OBJS="$OBJS ${src%.*}.o"
done

mkdir -p "$STELLA_ROOT/run/strad"
$FC -qopenmp -o "$STELLA_ROOT/run/strad/xstella6new.exe" $OBJS $MKL_LIBS 2>&1 | grep -v "remark #10448"

if [[ -f "$STELLA_ROOT/run/strad/xstella6new.exe" ]]; then
    ln -sf "$STELLA_ROOT/run/strad/xstella6new.exe" "$STELLA_ROOT/bin/xstella.exe"
    echo "SUCCESS: $STELLA_ROOT/bin/xstella.exe built"
else
    echo "ERROR: xstella6new.exe not found after link" >&2
    exit 1
fi

# --- Step 4: vladsf lineatom symlink ---------------------------------------
echo ""
echo "--- Step 4: vladsf lineatom.dat symlink ---"
cd "$STELLA_ROOT/run/vladsf"
if [[ ! -f lineatom.dat && -f lineatom_may07.dat ]]; then
    ln -sf lineatom_may07.dat lineatom.dat
fi

# --- Step 5: eve2 note ------------------------------------------------------
echo ""
echo "--- Step 5: eve2 (initial-model builder) ---"
if [[ ! -f "$STELLA_ROOT/run/eve/eve2.exe" ]]; then
    echo "NOTE: eve2.exe not present. Build it with the standard OpenStella"
    echo "      argmodelGF.sh flow if you need to generate .mod/.xni inputs."
else
    echo "eve2 already present"
fi

echo ""
echo "======================================"
echo "OpenStella_parallel build complete."
echo "Main executable : $STELLA_ROOT/bin/xstella.exe"
echo "Solver          : Intel MKL PARDISO (threaded)"
echo "======================================"
