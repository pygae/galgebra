# Release process runbook

This page is a step-by-step guide for cutting a galgebra release.
It captures lessons learned from the 0.5.1 incident (#517) and the 0.6.0 cycle.

## Policy

- Always cut a release candidate (RC) first.  The RC validates the full
  pipeline (PyPI upload, RTD build, changelog rendering) before the final tag.
- RC releases go to PyPI as pre-releases.  Zenodo skips pre-releases
  automatically; it archives only full releases.
- All PRs must be **squash-merged** (resumed since 0.6.0).

---

## Step 1 — Prepare the changelog

1. Ensure all PRs intended for the release are merged and have entries in
   `doc/changelog.rst`.

2. Add a `:release:` marker for the final version **only**.  Do **not** add
   `:release:` for the RC — `semantic_version` (used internally by the
   `releases` Sphinx extension) does not understand PEP 440 RC suffixes like
   `0.6.0rc1`, and the RTD build will fail.

   ```rst
   - :release:`0.6.0 <YYYY.MM.DD>`
   ```

3. Verify the `test/.nbval_sanitize.cfg` version regex covers all suffix
   patterns before tagging:

   ```
   regex: \{GAlgebra \}(\d+\.\d+(\.\d+)?(rc\d+|-dev)?)
   replace: {GAlgebra }
   ```

   If you add a new suffix pattern (e.g. `.post`), update this regex first.

---

## Step 2 — Bump the version for the RC

> **Critical:** the version string in `galgebra/_version.py` **must exactly
> match** the git tag you will push.  A mismatch (e.g. `_version.py` says
> `0.5.1` but the tag is `v0.5.1rc2`) caused the 0.5.1 Zenodo incident (#517)
> — PyPI accepted `0.5.1` from the RC tag and blocked the final release.

1. Create a branch:

   ```bash
   git checkout -b utensil/X.Y.ZrcN master
   ```

2. Edit `galgebra/_version.py`:

   ```python
   __version__ = 'X.Y.ZrcN'   # e.g. '0.6.0rc1'
   ```

3. Open a PR, get it merged (squash).

---

## Step 3 — Tag and create the RC release

```bash
git tag vX.Y.ZrcN
git push origin vX.Y.ZrcN

gh release create vX.Y.ZrcN \
    --prerelease \
    --title "vX.Y.ZrcN" \
    --notes "Release candidate for X.Y.Z.  See the [changelog](https://galgebra.readthedocs.io/en/latest/changelog.html) for details."
```

The CI `Create release and send to PyPI` job triggers automatically on the tag.

---

## Step 4 — Validate the RC

- **PyPI**: confirm the pre-release appears at https://pypi.org/project/galgebra/#history
- **RTD**: confirm the docs build passes and the changelog renders correctly
- **CI**: all Python versions green on the tag run
- **Install test**: `pip install galgebra==X.Y.ZrcN` and run a quick smoke test

---

## Step 5 — Bump the version for the final release

1. Create a branch:

   ```bash
   git checkout -b utensil/X.Y.Z master
   ```

2. Edit `galgebra/_version.py`:

   ```python
   __version__ = 'X.Y.Z'   # e.g. '0.6.0'
   ```

3. Add the `:release:` marker to `doc/changelog.rst` (see Step 1).

4. Open a PR, get it merged (squash).

---

## Step 6 — Tag and create the final release

```bash
git tag vX.Y.Z
git push origin vX.Y.Z

gh release create vX.Y.Z \
    --title "vX.Y.Z" \
    --notes "See the [changelog](https://galgebra.readthedocs.io/en/latest/changelog.html) for details."
```

---

## Step 7 — Post-release checks

| Check | How |
|-------|-----|
| PyPI | https://pypi.org/project/galgebra/#history |
| Zenodo | Verify a new record was created at https://zenodo.org/search?q=galgebra; if the webhook missed it, trigger manually from the Zenodo settings page for the record |
| Close milestone | `gh api repos/pygae/galgebra/milestones --jq '.[] | select(.title=="X.Y.Z") | .number'` then `gh api -X PATCH repos/pygae/galgebra/milestones/N -f state=closed` |
| README / docs | Open a follow-up PR to update any version references (badge, Prerequisites, install instructions) |
