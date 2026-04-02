# Bumping the Python version support

This page documents the process for changing the supported Python version
range and explains the policy behind it.

## Version policy

- **Drop** Python versions that have reached [end-of-life](https://devguide.python.org/versions/).
- **Add** the most recent minor release that is not the very latest (e.g. add
  3.12 when 3.13 is already out).  This keeps us current without chasing
  bleeding-edge releases.

The same "recent but not latest" rule applies to the SymPy version bump; see
[Bumping the SymPy dependency](bumping-sympy.md) for that workflow.

## Files to update

| File | What to change |
|------|----------------|
| `.github/workflows/ci.yml` | `python-version` matrix; runner OS; action versions; conditions using a specific version (e.g. `== '3.11'`) |
| `setup.py` | `python_requires`; `Programming Language :: Python :: X.Y` classifiers |
| `test_requirements.txt` | Add or update compatibility shims (e.g. `packaging` to replace `distutils`) |
| `README.md` | Prerequisites section Python version line and history note |

## CI action versions

When bumping the runner OS (e.g. `ubuntu-22.04` → `ubuntu-24.04`) also audit
the pinned action versions for OS compatibility:

- `actions/checkout` — bump to latest v4+
- `actions/setup-python` — bump to latest v5+
- `actions/cache` — bump to latest v4+
- `hidakatsuya/action-setup-diff-pdf` — check release notes; v1.4.0 added
  Ubuntu 24.04 support

## Python 3.12 compatibility audit

Before opening the bump PR, audit the codebase and notebooks for known Python
3.12 breaking changes:

**`distutils` removed**

```bash
grep -r "import distutils" .
```

Replace with `packaging` (add to `test_requirements.txt` if not present).

**Unrecognized escape sequences become `SyntaxWarning`**

In Python 3.12, `'\mathbf'`, `'\grad'`, etc. are `SyntaxWarning` (previously
`DeprecationWarning`).  nbval catches the stderr output and fails the cell.

```bash
python -W error::SyntaxWarning -c "import ast; ast.parse(open('file.py').read())"
```

For notebooks, search for unquoted backslash sequences in code cells and add
the `r` prefix: `'\mathbf{e}'` → `r'\mathbf{e}'`.

## Testing locally

```bash
uv venv --python 3.12 .venv312
uv pip install -r test_requirements.txt -e .
python -m flake8 -v
python -m pytest \
    -vv --durations=50 \
    --cov=galgebra \
    --nbval examples/ipython/ \
    --nbval examples/primer/ \
    test \
    --nbval-current-env \
    --nbval-sanitize-with test/.nbval_sanitize.cfg \
    -n 2 --dist loadscope
```

## Tracking issue pattern

Python bumps often uncover several independent blocking issues (compat
failures, notebook output changes).  Use the **tracking issue** pattern: one
tracker issue for the overall bump, sub-issues for each blocker, a PR per
sub-issue.  See [Bumping the SymPy dependency](bumping-sympy.md) for a
detailed description of this pattern.
