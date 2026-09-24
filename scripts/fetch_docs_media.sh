#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

# Fetch the figures and videos of the documentation into docs/media.
#
# They are kept in the repository scfpl-umich/COMPAS-media, and docs/media.lock records the
# commit of it that this version of the documentation uses. Run this before building the docs:
#
#     bash scripts/fetch_docs_media.sh
#     sphinx-build -b html docs docs/_build/html
#
# COMPAS_MEDIA_URL sets another address for the media repository, for example a local clone.

set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
lock="${root}/docs/media.lock"
dest="${root}/docs/media"
url="${COMPAS_MEDIA_URL:-https://github.com/scfpl-umich/COMPAS-media.git}"

commit=$(grep -v '^#' "${lock}" | tr -d '[:space:]')
if [ -z "${commit}" ]; then
    echo "No commit found in ${lock}" >&2
    exit 1
fi

if [ ! -d "${dest}/.git" ]; then
    git init -q "${dest}"
    git -C "${dest}" remote add origin "${url}"
else
    git -C "${dest}" remote set-url origin "${url}"
fi

if [ "$(git -C "${dest}" rev-parse -q --verify HEAD 2>/dev/null || true)" != "${commit}" ]; then
    git -C "${dest}" fetch -q --depth 1 origin "${commit}"
    git -C "${dest}" -c advice.detachedHead=false checkout -q --force FETCH_HEAD
fi

echo "docs/media is at ${commit}"
