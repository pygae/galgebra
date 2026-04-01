#!/usr/bin/env python3
"""Validate that a notebook re-execution introduced only cosmetic output changes.

Usage
-----
Compare two notebook files directly:

    python scripts/validate_nb_refresh.py OLD.ipynb NEW.ipynb

Use --git to extract the OLD version automatically from a git ref, comparing
it against the current working-tree file (NEW).  BEFORE_REF is any git
revision that existed before the notebook was re-executed:

    python scripts/validate_nb_refresh.py --git BEFORE_REF NOTEBOOK_PATH

Examples:

    # Compare the version at the previous commit (HEAD^) against the
    # current working copy -- typical use when reviewing a refresh PR:
    python scripts/validate_nb_refresh.py --git HEAD^ examples/ipython/gr_metrics.ipynb

    # Compare the version on master against the current branch:
    python scripts/validate_nb_refresh.py --git master examples/ipython/gr_metrics.ipynb

    # Compare an explicit saved copy against the refreshed file:
    git show HEAD^:examples/ipython/gr_metrics.ipynb > /tmp/old.ipynb
    python scripts/validate_nb_refresh.py /tmp/old.ipynb examples/ipython/gr_metrics.ipynb

Exits 0 if every output difference is accounted for by known cosmetic changes.
Exits 1 and prints a detailed report if any unexpected difference is found.

Intended use: run locally when reviewing a notebook-refresh PR to confirm the
re-execution changed nothing mathematically.  See doc/dev/bumping-sympy.rst for
the full workflow.

Known cosmetic changes handled
-------------------------------
1. LaTeX array column specs: ``\\begin{array}{ccc}`` <-> ``\\begin{array}{}``
   SymPy versions differ on whether they emit column-alignment characters.

2. LaTeX multiplication dot: ``\\cdot \\left(`` <-> ``\\left(``
   SymPy versions differ on whether an explicit ``\\cdot`` is printed when
   multiplying a radical by a parenthesized expression.

3. Plain-text whitespace: unicode matrix art uses spaces for alignment.
   Any whitespace-only difference is considered cosmetic.

4. Stream outputs (stdout/stderr): warnings from Python packages such as
   DeprecationWarnings from mpmath are environment-specific and ignored.
   Only ``display_data`` and ``execute_result`` outputs are compared.
"""

import json
import re
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Normalizers for text/latex outputs
# ---------------------------------------------------------------------------

def _norm_array_colspec(text: str) -> str:
    """``\\begin{array}{ccc}`` -> ``\\begin{array}{}`` (any column spec)."""
    return re.sub(r'\\begin\{array\}\{[^}]*\}', r'\\begin{array}{}', text)


def _norm_cdot(text: str) -> str:
    r"""``\cdot \left(`` -> ``\left(`` (explicit mult dot added/removed)."""
    return text.replace(r'\cdot \left(', r'\left(')


LATEX_NORMALIZERS = [
    _norm_array_colspec,
    _norm_cdot,
]


def normalize_latex(text: str) -> str:
    for fn in LATEX_NORMALIZERS:
        text = fn(text)
    return text


def normalize_plaintext(text: str) -> str:
    """Strip all whitespace — unicode matrix art uses spaces only for alignment."""
    return re.sub(r'\s+', '', text)


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

OUTPUT_TYPES_WITH_DATA = {'display_data', 'execute_result'}


def data_outputs(cell: dict) -> list[dict]:
    """Return only outputs with mathematical data, skipping stream outputs."""
    return [o for o in cell.get('outputs', []) if o.get('output_type') in OUTPUT_TYPES_WITH_DATA]


def compare_outputs(ci: int, old_outs: list, new_outs: list) -> list[str]:
    """Return failure messages for unexpected diffs between two output lists."""
    failures = []

    if len(old_outs) != len(new_outs):
        failures.append(
            f"  cell {ci}: data output count changed "
            f"({len(old_outs)} -> {len(new_outs)})"
        )
        return failures

    for oi, (oo, no) in enumerate(zip(old_outs, new_outs)):
        old_mimes = set(oo.get('data', {}).keys())
        new_mimes = set(no.get('data', {}).keys())

        if old_mimes != new_mimes:
            failures.append(
                f"  cell {ci} output {oi}: mime types changed {old_mimes} -> {new_mimes}"
            )
            continue

        for mime in old_mimes:
            ov = ''.join(oo['data'].get(mime, []))
            nv = ''.join(no['data'].get(mime, []))

            if ov == nv:
                continue

            # Apply normalization: text/plain can contain either unicode art
            # (normalize whitespace) or LaTeX strings (normalize latex patterns).
            # Apply all normalizers to text/plain since both formats appear there.
            if mime in ('text/latex', 'text/plain'):
                ov_n = normalize_plaintext(normalize_latex(ov))
                nv_n = normalize_plaintext(normalize_latex(nv))
            else:
                ov_n, nv_n = ov, nv

            if ov_n == nv_n:
                continue

            # Unexpected difference — report first differing line
            old_lines = ov.splitlines()
            new_lines = nv.splitlines()
            for lo, ln in zip(old_lines, new_lines):
                if lo != ln:
                    failures.append(
                        f"  cell {ci} output {oi} [{mime}]: unexpected diff\n"
                        f"    OLD: {lo[:120]!r}\n"
                        f"    NEW: {ln[:120]!r}"
                    )
                    break
            else:
                failures.append(
                    f"  cell {ci} output {oi} [{mime}]: unexpected diff "
                    f"(line count {len(old_lines)} -> {len(new_lines)})"
                )

    return failures


# ---------------------------------------------------------------------------
# Main validator
# ---------------------------------------------------------------------------

def validate(old_path: Path, new_path: Path) -> bool:
    with old_path.open() as f:
        old_nb = json.load(f)
    with new_path.open() as f:
        new_nb = json.load(f)

    old_code = [c for c in old_nb['cells'] if c['cell_type'] == 'code']
    new_code = [c for c in new_nb['cells'] if c['cell_type'] == 'code']

    if len(old_code) != len(new_code):
        print(f"FAIL {new_path.name}: code cell count changed "
              f"({len(old_code)} -> {len(new_code)})")
        return False

    all_failures = []
    for ci, (oc, nc) in enumerate(zip(old_code, new_code)):
        all_failures.extend(
            compare_outputs(ci, data_outputs(oc), data_outputs(nc))
        )

    if all_failures:
        print(f"FAIL {new_path.name}: {len(all_failures)} unexpected difference(s):")
        for msg in all_failures:
            print(msg)
        return False

    print(f"OK   {new_path.name}: all output diffs are cosmetic")
    return True


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = sys.argv[1:]

    if len(args) == 3 and args[0] == '--git':
        ref, nb_path_str = args[1], args[2]
        nb_path = Path(nb_path_str)
        result = subprocess.run(
            ['git', 'show', f'{ref}:{nb_path_str}'],
            capture_output=True,
        )
        if result.returncode != 0:
            print(f"Error: git show {ref}:{nb_path_str} failed:\n{result.stderr.decode()}")
            sys.exit(2)
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.ipynb', delete=False, mode='wb') as tmp:
            tmp.write(result.stdout)
            old_path = Path(tmp.name)
        ok = validate(old_path, nb_path)
        old_path.unlink()
        sys.exit(0 if ok else 1)

    if len(args) != 2:
        print(__doc__)
        sys.exit(2)

    old_path, new_path = Path(args[0]), Path(args[1])
    for p in (old_path, new_path):
        if not p.exists():
            print(f"Error: {p} does not exist")
            sys.exit(2)

    ok = validate(old_path, new_path)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
