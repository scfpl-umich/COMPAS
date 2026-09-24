# Contributing

Bug reports, fixes and new features are welcome.
[`CONTRIBUTING.md`](https://github.com/scfpl-umich/COMPAS/blob/main/CONTRIBUTING.md) is the full
guide. In short:

- Open a pull request against `main` from a branch in your fork. A maintainer reviews and approves
  it before it is merged, and the automated tests must pass.
- Run the test suite first, and extend it for what you add. The suite is described in
  [Testing](../user-guide/testing.md).
- A change that forces an edit under `exec/` changes the case interface. It must update every case
  and be recorded in `CHANGELOG.md`.
- Sign your commits with `git commit -s` (Developer Certificate of Origin). Contributions are
  licensed under the GPL, version 3 or later, like the rest of COMPAS (see [License](license.md)).

To report a bug, open an [issue](https://github.com/scfpl-umich/COMPAS/issues) with the details
that `CONTRIBUTING.md` lists.
