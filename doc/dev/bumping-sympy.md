# Bumping the SymPy dependency

This page documents the process for bumping the pinned SymPy version in CI
and explains the practices that have evolved around it.

## Version policy

Target the most recent **minor series** that is not the latest release.
For example, if 1.14 is the latest stable release, target 1.13.x (the
highest patch, e.g. 1.13.3).  This gives us a version that is current
enough to catch real compatibility issues but still widely used and stable.

The same policy applies to Python version bumps.

## Tracking issue pattern

A SymPy (or Python) bump typically uncovers several independent blocking
issues.  Use the **tracking issue** pattern to manage this:

1. Open a *tracker* issue describing the overall goal (e.g. "Bump SymPy
   to 1.13").  Keep it focused on the goal — do not list specific errors in
   the tracker body.
2. For each blocking problem found during testing, open a *sub-issue*
   describing only that problem.
3. Add a comment to the tracker referencing all sub-issues.
4. Open a PR for each sub-issue independently.  Sub-PRs can be reviewed
   and merged in any order.
5. Once all sub-PRs are merged, open the final bump PR that updates
   `test_requirements.txt` and closes the tracker.

This keeps each PR small and reviewable, and makes it easy to see at a
glance which problems are still blocking the bump.

## Testing locally

Use `uv` to create a clean venv with the target Python and SymPy versions:

```bash
uv venv --python 3.12 .venv312
uv pip install -r test_requirements.txt -e .
```

Run the full CI test command:

```bash
pytest \
    -vv --durations=50 \
    --cov=galgebra \
    --nbval examples/ipython/ \
    --nbval examples/primer/ \
    test \
    --nbval-current-env \
    --nbval-sanitize-with test/.nbval_sanitize.cfg \
    -n 2 --dist loadscope
```

CI also accepts a `PYTEST_K_FILTER` environment variable to pass `-k`.
See `.github/workflows/ci.yml` for the authoritative command.

## Common compatibility issues

### `distutils` removal (Python 3.12)

`distutils` was removed in Python 3.12.  Any `import distutils.*` must be
replaced with `packaging` (already a test dependency).

### Unrecognized escape sequences (Python 3.12)

Python 3.12 promotes `\m`, `\g`, `\i`, etc. in regular string literals from
`DeprecationWarning` to `SyntaxWarning`.  nbval catches the resulting stderr
output and fails the cell.  Fix: add the `r` prefix to affected string
literals (`r'\mathbf{e}'`).

### `ImmutableDenseMatrix` (SymPy 1.13)

Some SymPy 1.13 operations return `ImmutableDenseMatrix` where 1.12 returned
`MutableDenseMatrix`.  Any `isinstance(x, Matrix)` check that should cover
both must be changed to `isinstance(x, MatrixBase)` (imported from
`sympy.matrices`).

### LaTeX printing changes (SymPy 1.13)

SymPy 1.13 changed matrix LaTeX output in two ways:

- `\begin{array}{cccc}` → `\begin{array}{}` (column specs removed)
- Explicit `\cdot` multiplication dot added/removed in some expressions

These are cosmetic; the computed values are identical.  Fix: re-execute the
affected notebooks (see *Notebook refresh workflow* below).

### Unicode matrix art alignment (SymPy 1.13)

Plain-text pretty-printing of matrices shifts by one space in some rows.
Also cosmetic; handled by the notebook refresh.

## Notebook refresh workflow

When stored notebook outputs become stale due to a SymPy printing change:

1. Re-execute the affected notebooks with the target SymPy version:

   ```bash
   pip install jupyter nbconvert
   jupyter nbconvert --to notebook --execute --inplace \
       examples/ipython/gr_metrics.ipynb
   ```

2. Validate that every output change is purely cosmetic using
   `scripts/validate_nb_refresh.py`:

   ```bash
   # Compare the previous commit (HEAD^) against the current working copy:
   python scripts/validate_nb_refresh.py --git HEAD^ \
       examples/ipython/gr_metrics.ipynb
   ```

   The script normalizes all known cosmetic patterns (array column specs,
   `\cdot`, whitespace, stream warnings) and fails with a detailed diff if
   anything else changed.  Exit code 0 means all changes are cosmetic.

3. Verify the refreshed notebooks pass nbval:

   ```bash
   pytest --nbval examples/ipython/gr_metrics.ipynb \
       --nbval-current-env \
       --nbval-sanitize-with test/.nbval_sanitize.cfg
   ```

4. Include the `validate_nb_refresh.py` output in the PR description as
   evidence that no mathematical content changed.

### Adding new cosmetic normalizers

If a new SymPy version introduces a new cosmetic rendering difference, add a
normalizer function to `scripts/validate_nb_refresh.py`:

- For LaTeX diffs: add to `LATEX_NORMALIZERS` list.
- For plain-text diffs: update `normalize_plaintext`.
- Document the SymPy version and nature of the change in a comment next to
  the normalizer function.
