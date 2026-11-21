#!/usr/bin/env bash
# Build vepyr against a local datafusion-bio-functions checkout into a throwaway venv.
#
# Usage: build_venv.sh <functions-checkout> <venv-dir>
#
# The build runs in a dedicated vepyr worktree, $GATE_ROOT/vepyr-build, whose Cargo.toml
# is rewritten to a path dependency on <functions-checkout>. The wheel is installed
# non-editable, so the venv does not depend on that worktree afterwards. A build with a
# path dependency never qualifies a release.
#
# Environment (all optional):
#   GATE_ROOT       state directory        (default: $HOME/research/benchmarks/local-gate)
#   VEPYR_REPO      vepyr git checkout     (default: <this repository>/../vepyr)
#   VEPYR_REF       vepyr revision to build (default: origin/master)
#   GATE_TOOLCHAIN  Rust toolchain         (default: 1.91.1)
#   GATE_RUSTFLAGS  RUSTFLAGS              (default: -C target-cpu=native)
set -euo pipefail

usage() {
  echo "usage: build_venv.sh <functions-checkout> <venv-dir>" >&2
  exit 2
}
[ $# -eq 2 ] || usage
[ -f "$1/datafusion/bio-function-vep/Cargo.toml" ] || {
  echo "error: $1 is not a datafusion-bio-functions checkout" >&2
  usage
}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
REPO_ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
FUNCTIONS=$(cd "$1" && pwd -P)
mkdir -p "$(dirname "$2")"
VENV="$(cd "$(dirname "$2")" && pwd -P)/$(basename "$2")"

GATE_ROOT=${GATE_ROOT:-$HOME/research/benchmarks/local-gate}
VEPYR_REPO=${VEPYR_REPO:-$REPO_ROOT/../vepyr}
VEPYR_REF=${VEPYR_REF:-origin/master}
GATE_TOOLCHAIN=${GATE_TOOLCHAIN:-1.91.1}
GATE_RUSTFLAGS=${GATE_RUSTFLAGS:--C target-cpu=native}
BUILD="$GATE_ROOT/vepyr-build"

git -C "$VEPYR_REPO" rev-parse --git-dir > /dev/null 2>&1 || {
  echo "error: VEPYR_REPO=$VEPYR_REPO is not a git checkout" >&2
  exit 2
}

# The build worktree: create it once, otherwise drop the previous path dependency first.
if [ ! -e "$BUILD/.git" ]; then
  mkdir -p "$GATE_ROOT"
  GIT_LFS_SKIP_SMUDGE=1 git -C "$VEPYR_REPO" worktree add --detach "$BUILD" "$VEPYR_REF" > /dev/null
else
  git -C "$BUILD" checkout -q -- Cargo.toml Cargo.lock
  GIT_LFS_SKIP_SMUDGE=1 git -C "$BUILD" checkout -q --detach "$VEPYR_REF"
fi

count=$(grep -c '^datafusion-bio-function-vep = ' "$BUILD/Cargo.toml" || true)
[ "$count" -eq 1 ] || {
  echo "error: expected one datafusion-bio-function-vep dependency line in $BUILD/Cargo.toml, found $count" >&2
  exit 1
}
DEP="datafusion-bio-function-vep = { path = \"$FUNCTIONS/datafusion/bio-function-vep\", features = [\"cache-builder\"] }" \
  awk '/^datafusion-bio-function-vep = / { print ENVIRON["DEP"]; next } { print }' \
  "$BUILD/Cargo.toml" > "$BUILD/Cargo.toml.gate"
mv "$BUILD/Cargo.toml.gate" "$BUILD/Cargo.toml"

echo "building vepyr $(git -C "$BUILD" rev-parse --short HEAD) against functions $(git -C "$FUNCTIONS" rev-parse --short HEAD) into $VENV"
if ! (
  cd "$BUILD"
  env -u CONDA_PREFIX -u VIRTUAL_ENV \
    RUSTUP_TOOLCHAIN="$GATE_TOOLCHAIN" RUSTFLAGS="$GATE_RUSTFLAGS" CARGO_INCREMENTAL=0 \
    CARGO_TARGET_DIR="$GATE_ROOT/target/$(basename "$VENV")" UV_PROJECT_ENVIRONMENT="$VENV" \
    uv sync --frozen --no-cache --no-default-groups --no-editable --reinstall-package vepyr
) > "$VENV.build.log" 2>&1; then
  echo "error: build failed; see $VENV.build.log" >&2
  tail -20 "$VENV.build.log" >&2
  exit 1
fi

"$VENV/bin/python" -P - "$FUNCTIONS" "$BUILD" "$VENV" "$GATE_TOOLCHAIN" "$GATE_RUSTFLAGS" << 'EOF'
import glob
import hashlib
import json
import subprocess
import sys

functions, build, venv, toolchain, rustflags = sys.argv[1:]


def git(repo, *args):
    return subprocess.check_output(["git", "-C", repo, *args], text=True).strip()


extensions = glob.glob(f"{venv}/lib/python*/site-packages/vepyr/_core*.so")
if len(extensions) != 1:
    sys.exit(f"error: expected one vepyr extension in {venv}, found {len(extensions)}")
identity = {
    "functions_path": functions,
    "functions_sha": git(functions, "rev-parse", "HEAD"),
    "functions_dirty": bool(git(functions, "status", "--porcelain")),
    "vepyr_sha": git(build, "rev-parse", "HEAD"),
    "extension_sha256": hashlib.sha256(open(extensions[0], "rb").read()).hexdigest(),
    "toolchain": toolchain,
    "rustflags": rustflags,
}
json.dump(identity, open(f"{venv}/GATE_IDENTITY.json", "w"), indent=2)
print(json.dumps(identity, indent=2))
EOF
