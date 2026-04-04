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

5. SymPy 1.13 ``trigsimp(method='old')`` algebraic form differences
   (curvilinear coordinates example, ``examples/ipython/LaTeX.ipynb``):

   a. Pythagorean identity: ``sin²(η)+sinh²(ξ)`` <-> ``-cos²(η)+cosh²(ξ)``
   b. Power factoring: ``(X²+Y²)^{3/2}`` <-> ``X²√(X²+Y²)+Y²√(X²+Y²)``
   c. Spherical Laplacian/grad-wedge-B: collected ``\\frac{r²A+rB+C}{r²}``
      form <-> distributed ``A + B/r + C/r²`` form.
   d. Whitespace inside ``\\frac{...}`` arguments.
   e. Outer ``\\left(…\\right)`` wrapper before a basis blade.

   **Known remaining differences (require symbolic algebra to verify):**

   * Spherical curl ``e_r`` component: SymPy writes a parenthesised
     ``(A/tan + B - C/sin²)|sin|`` form; trigsimp produces a single
     fraction over ``tan|sin|``.  Both are mathematically equal but have
     fundamentally different structure.

   * Prolate-spheroidal divergence: SymPy collects all terms into one
     large fraction; trigsimp distributes into seven separate fractions.
     The two forms are mathematically equal and can be verified with
     SymPy but are not normalizable by simple string rewriting.
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


def _norm_sin2_sinh2_identity(text: str) -> str:
    r"""``{\\sin{ARG}}^{2} + {\\sinh{ARG2}}^{2}`` -> ``- {\\cos{ARG}}^{2} + {\\cosh{ARG2}}^{2}``.

    SymPy versions differ on which Pythagorean form they choose for the
    prolate-spheroidal metric coefficient sin²(η) + sinh²(ξ) = cosh²(ξ) - cos²(η).
    """
    pattern = (
        r'\{\\sin\{(\\left \(.*?\\right \))\}\}\^\{2\}'
        r' \+ '
        r'\{\\sinh\{(\\left \(.*?\\right \))\}\}\^\{2\}'
    )
    return re.sub(pattern, r'- {\\cos{\1}}^{2} + {\\cosh{\2}}^{2}', text)


def _norm_sum_sqrt_to_power32(text: str) -> str:
    r"""``\\left(X^{2} + Y^{2}\\right)^{\\frac{3}{2}}`` -> ``X^{2}\\sqrt{...} + Y^{2}\\sqrt{...}``.

    SymPy versions differ on whether (u²+v²)^{3/2} is written as a single power
    or factored as u²·√(u²+v²) + v²·√(u²+v²).
    """
    pattern = (
        r'\\left\(([a-z])\^\{2\} \+ ([a-z])\^\{2\}\\right\)'
        r'\^\{\\frac\{3\}\{2\}\}'
    )
    replacement = r'\1^{2} \\sqrt{\1^{2} + \2^{2}} + \2^{2} \\sqrt{\1^{2} + \2^{2}}'
    return re.sub(pattern, replacement, text)


# ---------------------------------------------------------------------------
# Brace-counting helpers for structural LaTeX normalizers
# ---------------------------------------------------------------------------

def _find_matching_brace(text: str, pos: int) -> int:
    """Return the index of the ``}`` matching the ``{`` at *pos*, or -1."""
    depth = 0
    for i in range(pos, len(text)):
        c = text[i]
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return i
    return -1


def _split_signed_terms(text: str) -> list:
    r"""Split *text* at top-level ``+``/``-`` separators.

    Returns a list of ``(sign, term)`` tuples where *sign* is ``'+'`` or
    ``'-'`` and *term* is the stripped content.  A leading ``'- '`` is
    treated as the sign of the first term.
    """
    text = text.strip()
    result = []
    depth = 0
    sign = '+'
    start = 0

    if text.startswith('- '):
        sign = '-'
        start = 2

    i = start
    while i < len(text):
        c = text[i]
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
        elif (depth == 0 and c in '+-'
              and i > 0 and text[i - 1] == ' '
              and i + 1 < len(text) and text[i + 1] == ' '):
            term = text[start:i].strip()
            if term:
                result.append((sign, term))
            sign = c
            start = i + 2
            i += 2
            continue
        i += 1

    term = text[start:].strip()
    if term:
        result.append((sign, term))
    return result


def _divide_term_by_r2(sign: str, term: str) -> tuple:
    r"""Return ``(sign, term/r^{2})`` in LaTeX, possibly adjusting *sign*."""
    # r^{2} STUFF  →  STUFF
    if term.startswith('r^{2} '):
        return sign, term[6:]

    # N r STUFF  (integer coefficient × r × rest)  →  \frac{N STUFF}{r}
    m = re.match(r'^(\d+) r (.+)', term)
    if m:
        return sign, r'\frac{' + m.group(1) + ' ' + m.group(2) + r'}{r}'

    # \frac{NUM}{DEN}  →  \frac{NUM}{r^{2} DEN}
    if term.startswith(r'\frac{'):
        num_end = _find_matching_brace(term, 5)
        if num_end != -1 and num_end + 1 < len(term) and term[num_end + 1] == '{':
            denom_end = _find_matching_brace(term, num_end + 1)
            if denom_end != -1:
                num = term[6:num_end]
                den = term[num_end + 2:denom_end]
                return sign, r'\frac{' + num + r'}{r^{2} ' + den + '}'

    # plain STUFF  →  \frac{STUFF}{r^{2}}
    return sign, r'\frac{' + term + r'}{r^{2}}'


def _norm_distribute_r2_denominator(text: str) -> str:
    r"""Distribute ``\\frac{r^{2} A + 2r B + C}{r^{2}}`` into expanded form.

    Converts the collected-fraction form that older SymPy produces to the
    term-by-term form produced by newer SymPy (or vice-versa).  Only fires
    when the denominator is exactly ``{r^{2}}`` and the numerator starts
    with ``r^{2} ``.
    """
    result: list = []
    i = 0
    while i < len(text):
        if text[i:i + 6] != r'\frac{':
            result.append(text[i])
            i += 1
            continue

        open_brace = i + 5          # index of the opening {
        num_end = _find_matching_brace(text, open_brace)
        if num_end == -1:
            result.append(text[i]); i += 1; continue

        # Denominator must be exactly {r^{2}}  (7 chars)
        if text[num_end + 1:num_end + 8] != '{r^{2}}':
            result.append(text[i]); i += 1; continue

        numerator = text[open_brace + 1:num_end]

        # Only distribute when numerator opens with r^{2}
        if not numerator.startswith('r^{2} '):
            result.append(text[i]); i += 1; continue

        terms = _split_signed_terms(numerator)
        parts: list = []
        for j, (sign, term) in enumerate(terms):
            _, new_term = _divide_term_by_r2(sign, term)
            if j == 0:
                parts.append(('- ' if sign == '-' else '') + new_term)
            else:
                parts.append(' ' + sign + ' ' + new_term)
        result.append(''.join(parts))
        i = num_end + 8     # skip past }{r^{2}}
    return ''.join(result)


def _norm_strip_outer_parens_before_basis(text: str) -> str:
    r"""Remove a lone ``\\left ( X\\right ) \\boldsymbol`` wrapper.

    SymPy sometimes wraps a distributed sum in ``\\left ( ... \\right )``
    before a basis blade; other versions omit this wrapping.
    The inner ``\\left (`` / ``\\right )`` pairs inside function arguments
    (e.g. ``\\tan{\\left (\\theta \\right )}``) are always followed by ``}``
    rather than `` \\boldsymbol``, so the non-greedy match stops at the
    correct level.

    Assumption: galgebra's LaTeX output never contains a bare nested
    ``\\left ( ... \\right )`` group that is itself immediately followed by
    `` \\boldsymbol``.  If it did, the non-greedy ``.*?`` would greedily stop
    at the *inner* ``\\right )`` and produce a wrong result.  This holds for
    all current galgebra output formats; revisit if new expression types are
    added that wrap sub-expressions in bare parentheses before a basis blade.
    """
    return re.sub(
        r'\\left \( (.*?)\\right \) \\boldsymbol',
        r'\1 \\boldsymbol',
        text,
        flags=re.DOTALL,
    )


def _norm_collapse_spaces(text: str) -> str:
    r"""Normalize insignificant whitespace in LaTeX math.

    SymPy emits trailing spaces inside ``\\frac{...}`` arguments and between
    terms (e.g. ``\\frac{2 f }{r}`` vs ``\\frac{2 f}{r}``).  In LaTeX math,
    a space before ``}`` is invisible.  This normalizer:

    * Collapses runs of two or more spaces to one space.
    * Strips spaces immediately before ``}``.
    """
    text = re.sub(r'  +', ' ', text)
    text = re.sub(r' +\}', '}', text)
    return text


LATEX_NORMALIZERS = [
    _norm_array_colspec,
    _norm_cdot,
    # SymPy 1.13 (trigsimp method='old') uses different algebraic forms for
    # some curvilinear-coordinate expressions.  The normalizers below bring
    # both forms to a common representation so the validator can confirm the
    # changes are cosmetic.
    _norm_sin2_sinh2_identity,      # sin²+sinh²  ↔  -cos²+cosh²  (prolate spheroidal)
    _norm_sum_sqrt_to_power32,      # (X²+Y²)^{3/2}  ↔  X²√(…)+Y²√(…)  (paraboloidal)
    _norm_distribute_r2_denominator,  # \frac{r²A+rB+C}{r²}  ↔  A+B/r+C/r²  (spherical)
    _norm_strip_outer_parens_before_basis,  # \left( X\right) basis  ↔  X basis
    _norm_collapse_spaces,          # trailing/extra spaces in \frac args
]


def normalize_latex(text: str) -> str:
    for fn in LATEX_NORMALIZERS:
        text = fn(text)
    return text


_BOX_DRAWING = re.compile(r'[\u2500-\u257F]')


def normalize_plaintext(text: str) -> str:
    """Normalize whitespace in unicode matrix-art lines only.

    Lines containing box-drawing characters (U+2500–U+257F) may shift
    alignment between SymPy versions — treat any whitespace-only diff there
    as cosmetic.  Lines without box-drawing chars are left verbatim so that
    a real content change (e.g. ``- x`` → ``-x``) is still caught.
    """
    lines = []
    for line in text.splitlines():
        if _BOX_DRAWING.search(line):
            lines.append(re.sub(r'\s+', '', line))
        else:
            lines.append(line)
    return '\n'.join(lines)


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
