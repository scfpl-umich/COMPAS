# Contributing to COMPAS

Bug reports, fixes and new features are welcome. This page says how changes get into the code.

## How changes land

- `main` is the only long-lived branch. Releases are tags on `main` (`v1.0.0`, `v1.1.0`, ...).
- Work on a branch in your fork and open a pull request against `main`.
- Every pull request needs an approving review from a maintainer before it can be merged, and the
  test suite must pass. Maintainers are requested automatically through `.github/CODEOWNERS`.
- Keep a pull request to one change. Pull requests are squash-merged into one commit on `main`,
  whose message is the pull request title followed by the commit messages and their sign-offs.

## Before you open a pull request

1. Build and run the test suite:

       bash scripts/test_cases.sh

   It builds every case in `exec/_Tests` and runs it and its variants briefly.
2. A new run-time option needs a line in the `variants.txt` of a test case for each model it
   applies to. A new model or compile-time option needs a case in `exec/_Tests/`.
3. Document new options in `docs/user-guide/inputs.md` (run time) or
   `docs/user-guide/compile-options.md` (compile time).
4. Start every new source file with the license header below.

## Changes that affect the cases

A change that forces an edit under `exec/` is a change to the case interface. That includes the
signature of a function in `ProblemICBC.H`, a member of `Parm`, a hook macro, the output helper in
`include/Tools/UserOutput.H`, or the name or meaning of an inputs key. Such a pull request must:

- update every case under `exec/` and `exec/_Tests/` in the same pull request, and
- add an entry under "Case interface" in `CHANGELOG.md` with a short before-and-after example, so
  users can update their own cases.

## Sign off your commits

COMPAS uses the Developer Certificate of Origin 1.1 (<https://developercertificate.org>). By adding a
`Signed-off-by` line to a commit you certify it: in short, that you wrote the change or otherwise
have the right to submit it under the license indicated in the file, and that you understand the
contribution and your sign-off, including your name and email, are public and kept permanently.
Add the line with

    git commit -s

Contributions are licensed under the GNU General Public License, version 3 or later, with the
additional permission under section 7 given in `COPYRIGHT`, like the rest of COMPAS. By
contributing you grant that permission for your contribution.

## License header for new files

```cpp
// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
```

In Python and shell scripts, `GNUmakefile` and `Make.package` files, use the same lines with `#`
in place of `//`, and put them after the `#!` line if the file has one. A file that contains code
derived from AMReX adds one more line after the header:

```cpp
// Portions derived from AMReX, BSD-3-Clause. See THIRD_PARTY_NOTICES.md.
```

## Documentation figures and videos

The figures and videos of the documentation are kept in
[scfpl-umich/COMPAS-media](https://github.com/scfpl-umich/COMPAS-media), so that the history of
this repository stays small. `docs/media.lock` records the commit of that repository the
documentation uses, and `scripts/fetch_docs_media.sh` fetches it into `docs/media` before a build:

    bash scripts/fetch_docs_media.sh
    sphinx-build -b html docs docs/_build/html

To add or replace a figure or video, first open a pull request in COMPAS-media. After it is
merged, open a pull request here that sets `docs/media.lock` to the new commit and updates the pages
that use the file. Pages refer to figures as `media/figures/<name>` and to videos as
`_static/<name>`.

## Reporting a bug

Open an issue with the COMPAS version or commit, the AMReX version, the compiler and MPI library,
the case's `GNUmakefile` and `inputs` file, and the output that shows the problem.
