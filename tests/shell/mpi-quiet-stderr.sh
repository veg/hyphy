#!/bin/bash
# Regression test for veg/hyphy#2035 (part 2): HYPHYMPI must not write ANSI
# status-line progress ("ESC[2K ... Working on site pattern") to a non-terminal
# stderr. mpirun always interposes a pipe on rank 0's stderr, so the binary can
# not rely on isatty(); it uses TERM (and an explicit QUIET switch) instead.
#
# usage: tests/shell/mpi-quiet-stderr.sh [path/to/HYPHYMPI] [path/to/mpirun]
#
# mpirun (OpenMPI 5/prrte) forwards rank 0's stderr, and depending on the
# launcher version the forwarded progress can surface on mpirun's stdout or
# stderr, so both streams are scanned for ESC[2K.
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
HYPHYMPI=${1:-$ROOT/HYPHYMPI}
MPIRUN=${2:-mpirun}
DATA=$ROOT/tests/hbltests/libv3/data/CD2.nex
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
fail=0
ESC=$'\033'
nesc() { cat "$WORK/$1.out" "$WORK/$1.err" | grep -c -F "${ESC}[2K"; }

run() { # run <label> <env-prefix...> -- <extra hyphy args...>
  local label=$1; shift
  local envs=(); while [ "$1" != "--" ]; do envs+=("$1"); shift; done; shift
  env "${envs[@]}" "$MPIRUN" --oversubscribe --allow-run-as-root -np 2 "$HYPHYMPI" \
      LIBPATH="$ROOT/res" "$@" fel --alignment "$DATA" --code Universal --branches All --srv No \
      --output "$WORK/$label.json" >"$WORK/$label.out" 2>"$WORK/$label.err" </dev/null
  rc=${PIPESTATUS[0]}; echo "[$label] exit=$rc stderr_bytes=$(wc -c <"$WORK/$label.err") esc2K=$(nesc $label)"
}


# 1. batch-like: TERM unset, stderr is a file => no progress lines
run noterm -u TERM --
[ "$(nesc noterm)" -eq 0 ] || { echo "FAIL: TERM unset leaked ANSI progress to stderr"; fail=1; }

# 2. TERM=dumb => quiet
run dumb TERM=dumb --
[ "$(nesc dumb)" -eq 0 ] || { echo "FAIL: TERM=dumb leaked ANSI progress to stderr"; fail=1; }

# 3. existing explicit switch (SUPPRESS_STATUS_UPDATES) wins even for an 'interactive' TERM
run quiet TERM=xterm -- ENV="SUPPRESS_STATUS_UPDATES=1;"
[ "$(nesc quiet)" -eq 0 ] || { echo "FAIL: ENV=SUPPRESS_STATUS_UPDATES=1 leaked ANSI progress to stderr"; fail=1; }

# 4. interactive-ish (TERM=xterm; mpirun still pipes stderr) keeps progress
run xterm TERM=xterm --
[ "$(nesc xterm)" -gt 0 ] || { echo "FAIL: TERM=xterm lost interactive progress display"; fail=1; }

[ $fail -eq 0 ] && echo "PASS" || echo "FAILED"
exit $fail
