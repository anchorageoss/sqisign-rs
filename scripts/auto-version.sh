#!/usr/bin/env bash
set -euo pipefail

# Compute a semantic version for a build/release.
#
# The MAJOR.MINOR pair is the authoritative release line: the script reads the
# full SemVer `[workspace.package] version` in Cargo.toml (e.g. "0.3.0") and
# keeps only its MAJOR.MINOR prefix ("0.3"). The PATCH component is derived from
# git distance, so on a clean `main` it is the commit height -- producing tags
# like `v0.3.<distance>` that stay on the crates.io `0.3.x` line. To start a new
# minor/major line, bump the version in Cargo.toml and the tags follow.

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
# the truncated history). On CI, deepen to full history and fail if we can't;
# locally, warn loudly so a bogus version is obvious rather than silent.
# Workflows that invoke this script should check out with fetch-depth: 0; this
# guard is belt-and-suspenders for when they don't.
if [ "$(git rev-parse --is-shallow-repository 2>/dev/null)" = "true" ]; then
  if [ "${GITHUB_ACTIONS:-}" = "true" ]; then
    echo "Repository is shallow; fetching full history for an accurate version..." >&2
    git fetch --unshallow 2>/dev/null || git fetch --depth=2147483647 2>/dev/null || true
    # Don't silently proceed on a still-shallow repo: a truncated commit height
    # would mint a bogus version/tag.
    if [ "$(git rev-parse --is-shallow-repository 2>/dev/null)" = "true" ]; then
      echo "ERROR: repository is still shallow after fetch; refusing to compute a version from truncated history." >&2
      exit 1
    fi
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
  # printf, not echo: a branch name starting with `-` (e.g. `-n`) would be
  # treated as an option by echo and corrupt the metadata.
  # shellcheck disable=SC2001
  BRANCH_META="$(printf '%s' "$RAW_BRANCH" | sed 's/[^a-zA-Z0-9-]/-/g')-"
else
  BRANCH_META=
fi

# Treat the build as a default-branch build only when we know the ref is
# actually the canonical default branch. PR events always represent a
# non-default ref by definition, so even a fork PR opened from a branch
# literally named `main`/`master` must fall through to the merge-base path
# below. Both `pull_request` and `pull_request_target` set GITHUB_HEAD_REF, so a
# non-empty GITHUB_HEAD_REF is the reliable "this is a PR" signal; the explicit
# event-name checks are belt-and-suspenders.
if { [ "$RAW_BRANCH" = "main" ] || [ "$RAW_BRANCH" = "master" ]; } \
   && [ -z "${GITHUB_HEAD_REF:-}" ] \
   && [ "${GITHUB_EVENT_NAME:-}" != "pull_request" ] \
   && [ "${GITHUB_EVENT_NAME:-}" != "pull_request_target" ]; then
  HEIGHT=$(git rev-list --count HEAD)
  echo "$BASE_VERSION.$HEIGHT+${BRANCH_META}$SHORT_HASH"
  exit 0
fi

# Find the remote that points at the canonical upstream repo. Prefer
# $GITHUB_REPOSITORY (set in Actions) so forks/renames work without editing this
# script; fall back to the known upstream slug, then to "origin". Override with
# AUTO_VERSION_REMOTE if your remote layout differs.
EXPECTED_REPO="${AUTO_VERSION_REPO:-${GITHUB_REPOSITORY:-anchorageoss/sqisign-rs}}"
REMOTE="${AUTO_VERSION_REMOTE:-}"
if [ -z "$REMOTE" ]; then
  REMOTE=$(git remote -v | awk -v repo="$EXPECTED_REPO" '$0 ~ "[[:space:]]\\(fetch\\)" && index($0, repo) {print $1; exit}')
fi
if [ -z "$REMOTE" ]; then
  REMOTE="origin"
fi

# Resolve the default branch. Prefer main, then master. If neither
# remote-tracking ref is present (e.g. a CI checkout that fetched only the
# current ref), try to fetch main before giving up, so the merge-base below
# doesn't fail cryptically under `set -e`.
if git rev-parse --verify "$REMOTE/main" > /dev/null 2>&1; then
  DEFAULT_BRANCH="main"
elif git rev-parse --verify "$REMOTE/master" > /dev/null 2>&1; then
  DEFAULT_BRANCH="master"
else
  # Try to fetch both candidate default branches (separately, so a missing
  # `master` ref doesn't abort the `main` fetch), then re-resolve.
  git fetch "$REMOTE" main > /dev/null 2>&1 || true
  git fetch "$REMOTE" master > /dev/null 2>&1 || true
  if git rev-parse --verify "$REMOTE/main" > /dev/null 2>&1; then
    DEFAULT_BRANCH="main"
  elif git rev-parse --verify "$REMOTE/master" > /dev/null 2>&1; then
    DEFAULT_BRANCH="master"
  else
    echo "ERROR: cannot resolve $REMOTE/main or $REMOTE/master to compute a merge base." >&2
    exit 1
  fi
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
