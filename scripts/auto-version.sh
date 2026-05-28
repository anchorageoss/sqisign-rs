#!/usr/bin/env bash
set -euo pipefail

# Compute a semantic version for a build/release.
#
# The MAJOR.MINOR pair is the authoritative release line, read from the
# workspace `[workspace.package] version` in Cargo.toml (e.g. "0.3"). The PATCH
# component is derived from git distance, so on a clean `main` it is the commit
# height -- producing tags like `v0.3.<distance>` that stay on the crates.io
# `0.3.x` line. To start a new minor/major line, bump the version in Cargo.toml
# and the tags follow automatically.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${SCRIPT_DIR}/../Cargo.toml"

# Read MAJOR.MINOR from the [workspace.package] version field. Only the version
# inside the [workspace.package] table is considered, so dependency `version =`
# lines elsewhere in the manifest can't be picked up by mistake.
BASE_VERSION="$(awk '
  /^\[workspace\.package\]/ { inblock = 1; next }
  /^\[/                     { inblock = 0 }
  inblock && /^[[:space:]]*version[[:space:]]*=/ {
    if (match($0, /"[^"]+"/)) {
      v = substr($0, RSTART + 1, RLENGTH - 2)
      split(v, a, ".")
      print a[1] "." a[2]
      exit
    }
  }
' "$MANIFEST")"

if [ -z "$BASE_VERSION" ]; then
  echo "ERROR: could not read [workspace.package] version from $MANIFEST" >&2
  exit 1
fi

# A shallow clone makes the commit-distance counts below wrong (git only sees
# the truncated history). On CI, deepen to full history; locally, warn loudly so
# a bogus version is obvious rather than silent. Callers in this repo check out
# with fetch-depth: 0, so this is a belt-and-suspenders guard.
if [ "$(git rev-parse --is-shallow-repository 2>/dev/null)" = "true" ]; then
  if [ "${GITHUB_ACTIONS:-}" = "true" ]; then
    echo "Repository is shallow; fetching full history for an accurate version..." >&2
    git fetch --unshallow >&2 2>/dev/null || git fetch --depth=2147483647 >&2 || true
  else
    echo "WARNING: shallow clone detected -- the computed version will be wrong." >&2
    echo "         Run 'git fetch --unshallow' for an accurate commit distance." >&2
  fi
fi

# Resolve the raw branch name before sanitization so comparisons against
# "main"/"master" can't be fooled by branches like "main." that collapse
# to "main-" after character sanitization.
#
# Precedence:
#   1. GITHUB_HEAD_REF -- set on pull_request events
#   2. GITHUB_REF_NAME -- set on push/workflow_dispatch events
#   3. git symbolic-ref --short HEAD -- local developer runs
RAW_BRANCH="${GITHUB_HEAD_REF:-}"
if [ -z "$RAW_BRANCH" ]; then
  RAW_BRANCH="${GITHUB_REF_NAME:-}"
fi
if [ -z "$RAW_BRANCH" ]; then
  if git symbolic-ref --short HEAD > /dev/null 2>&1; then
    RAW_BRANCH="$(git symbolic-ref --short HEAD)"
  fi
fi

SHORT_HASH=$(git rev-parse --short=12 HEAD)

# Sanitized branch for semver build metadata (allowed: [0-9A-Za-z-]).
# Trailing "-" is a separator before SHORT_HASH; omitted when branch is empty.
if [ -n "$RAW_BRANCH" ]; then
  # shellcheck disable=SC2001
  BRANCH_META="$(echo "$RAW_BRANCH" | sed 's/[^a-zA-Z0-9-]/-/g')-"
else
  BRANCH_META=
fi

# Treat the build as a default-branch build only when we know the ref is
# actually the canonical default branch. `pull_request` events always
# represent a non-default ref by definition, so even a fork PR opened from
# a branch literally named `main`/`master` falls through to the merge-base
# path below.
if { [ "$RAW_BRANCH" = "main" ] || [ "$RAW_BRANCH" = "master" ]; } \
   && [ "${GITHUB_EVENT_NAME:-}" != "pull_request" ]; then
  HEIGHT=$(git rev-list --count HEAD)
  echo "$BASE_VERSION.$HEIGHT+${BRANCH_META}$SHORT_HASH"
  exit 0
fi

REMOTE=$(git remote -v | awk '/[[:space:]]\(fetch\)/ && /anchorageoss\/sqisign-rs/ {print $1; exit}')
if [ -z "$REMOTE" ]; then
  REMOTE="origin"
fi

DEFAULT_BRANCH="main"
if ! git rev-parse --verify "$REMOTE/$DEFAULT_BRANCH" > /dev/null 2>&1; then
  DEFAULT_BRANCH="master"
fi

MERGE_BASE=$(git merge-base "$REMOTE/$DEFAULT_BRANCH" HEAD)
if [ "$MERGE_BASE" = "$(git rev-parse "$REMOTE/$DEFAULT_BRANCH")" ] && [ "${GITHUB_ACTIONS:-}" = "true" ]; then
  # On CI, the remote-tracking ref may be stale (shallow clone) -- fetch to get the real merge base.
  # Skipped on local builds; run `git fetch` manually if the version number looks wrong.
  echo "Fetching $REMOTE..." >&2
  git fetch "$REMOTE" >&2
  MERGE_BASE=$(git merge-base "$REMOTE/$DEFAULT_BRANCH" HEAD)
fi
# Patch = the main commit height at the merge-base: a stable base number shared
# with the main line, not the branch's own distance. The commit count since the
# merge-base (MERGE_DIFF) is carried in build metadata so distinct branch states
# stay distinguishable; only main mints release tags, so branch patches need not
# be globally unique.
MERGE_HEIGHT=$(git rev-list --count "$MERGE_BASE")
HEIGHT=$(git rev-list --count HEAD)
MERGE_DIFF=$((HEIGHT - MERGE_HEIGHT))
echo "$BASE_VERSION.$MERGE_HEIGHT+${BRANCH_META}${MERGE_DIFF}-$SHORT_HASH"
