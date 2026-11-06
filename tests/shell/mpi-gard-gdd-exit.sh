#!/bin/bash
# Regression test for veg/hyphy#2016: GARD under HYPHYMPI wrote a complete JSON
# and then never exited (rank 0 spinning, workers waiting for the shutdown
# message). The master was stuck in GARD.bf, not in MPI: the multi-breakpoint
# GA feasibility guard was a no-op (a 'gard.'-qualified name inside
# 'namespace gard' resolves to gard.gard.*, i.e. 0), so the GA started with a
# number of breakpoints that cannot fit, and gard.GA.initializeModels
# rejection-sampled forever.
#
# The bundled alignment (60 sequences x 300 nt => minimum partition size 117)
# admits one breakpoint (2*117 <= 300) but not two (3*117 > 300), which is
# exactly the hanging regime (the GA is entered after the single-breakpoint
# scan, then cannot place 2 breakpoints). Both the reported --rv GDD case and
# the --rv None control are run; each must exit 0 within TIMEOUT seconds and
# leave a parseable JSON. Runtime with the fix: ~30 s at NP=8.
#
# usage: tests/shell/mpi-gard-gdd-exit.sh [path/to/HYPHYMPI] [path/to/mpirun]
# env:   NP (ranks, default 4), TIMEOUT (seconds per run, default 300),
#        RES (LIBPATH, default <repo>/res)
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
HYPHYMPI=${1:-$ROOT/HYPHYMPI}
MPIRUN=${2:-mpirun}
NP=${NP:-4}
TIMEOUT=${TIMEOUT:-300}
RES=${RES:-$ROOT/res}
DATA=$ROOT/tests/hbltests/libv3/data/gard-sim60x300.fas
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
fail=0

run() { # run <label> <extra GARD args...>
  local label=$1; shift
  # GARD writes <alignment>.best-gard next to the input, so work on a copy
  cp "$DATA" "$WORK/$label.fas"
  local t0=$SECONDS
  timeout -k 10 "$TIMEOUT" "$MPIRUN" --oversubscribe --allow-run-as-root -np "$NP" "$HYPHYMPI" \
      LIBPATH="$RES" GARD --type nucleotide --alignment "$WORK/$label.fas" --mode Faster "$@" \
      --output "$WORK/$label.json" >"$WORK/$label.out" 2>"$WORK/$label.err" </dev/null
  local rc=$?
  echo "[$label] exit=$rc wall=$((SECONDS - t0))s"
  if [ $rc -eq 124 ] || [ $rc -eq 137 ]; then
    echo "FAIL: [$label] did not exit within ${TIMEOUT}s (hang, veg/hyphy#2016)"; fail=1
  elif [ $rc -ne 0 ]; then
    echo "FAIL: [$label] exited with $rc"; tail -5 "$WORK/$label.err" "$WORK/$label.out"; fail=1
  fi
  if ! python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert "breakpointData" in d' "$WORK/$label.json" 2>/dev/null; then
    echo "FAIL: [$label] output JSON missing or unparseable"; fail=1
  fi
}

run gdd  --rv GDD --rate-classes 3
run none --rv None

[ $fail -eq 0 ] && echo "PASS" || echo "FAILED"
exit $fail
