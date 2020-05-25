r"""
ANSI Enhanced Text Printing, Text Printer and LaTeX Printer for all Geometric Algebra classes
"""

import copy
import os
import sys
import io
import builtins
import functools
import inspect
import re
import shutil
from collections import ChainMap

from sympy import Matrix, Basic, S, Symbol, Function, Derivative, Pow
from sympy.printing.str import StrPrinter
from sympy.printing.conventions import split_super_sub
from sympy.printing.latex import LatexPrinter, accepted_latex_functions
from sympy.core.function import _coeff_isneg
from sympy.core.operations import AssocOp
from sympy import init_printing

try:
    from IPython.display import display, Latex, Math, display_latex
except ImportError:
    pass
try:
    from sympy.interactive import printing
except ImportError:
    pass

from inspect import getouterframes, currentframe

from ._utils import parser as _parser

ZERO_STR = ' 0 '

Format_cnt = 0

ip_cmds = r"""
$\DeclareMathOperator{\Tr}{Tr}
\DeclareMathOperator{\Adj}{Adj}
\newcommand{\bfrac}[2]{\displaystyle\frac{#1}{#2}}
\newcommand{\lp}{\left (}
\newcommand{\rp}{\right )}
\newcommand{\paren}[1]{\lp {#1} \rp}
\newcommand{\half}{\frac{1}{2}}
\newcommand{\llt}{\left <}
\newcommand{\rgt}{\right >}
\newcommand{\abs}[1]{\left |{#1}\right | }
\newcommand{\pdiff}[2]{\bfrac{\partial {#1}}{\partial {#2}}}
\newcommand{\npdiff}[3]{\bfrac{\partial^{#3} {#1}}{\partial {#2}^{#3}}}
\newcommand{\lbrc}{\left \{}
\newcommand{\rbrc}{\right \}}
\newcommand{\W}{\wedge}
\newcommand{\prm}[1]{{#1}'}
\newcommand{\ddt}[1]{\bfrac{d{#1}}{dt}}
\newcommand{\R}{\dagger}
\newcommand{\deriv}[3]{\bfrac{d^{#3}#1}{d{#2}^{#3}}}
\newcommand{\grade}[1]{\left < {#1} \right >}
\newcommand{\f}[2]{{#1}\lp {#2} \rp}
\newcommand{\eval}[2]{\left . {#1} \right |_{#2}}$
"""

SYS_CMD = {'linux2': {'rm': 'rm', 'evince': 'evince', 'null': ' > /dev/null', '&': '&'},
           'linux': {'rm': 'rm', 'evince': 'evince', 'null': ' > /dev/null', '&': '&'},
           'win32': {'rm': 'del', 'evince': 'start', 'null': ' > NUL', '&': ''},
           'darwin': {'rm': 'rm', 'evince': 'open', 'null': ' > /dev/null', '&': '&'}}


def isinteractive():  #Is ipython running
    """
    We will assume that if ipython is running then jupyter notebook is
    running.
    """
    try:
        __IPYTHON__
        return True
    except NameError:
        return False


def ostr(obj, dict_mode=False, indent=True):
    return GaPrinter(dict(dict_mode=dict_mode)).doprint(obj)


def find_functions(expr):
    f_lst = []
    for f in list(expr.atoms(Function)):
        if str(f) not in GaPrinter.function_names:
            f_lst.append(f)
    f_lst += list(expr.atoms(Derivative))
    return f_lst


def coef_simplify(expr):
    """
    fcts = find_functions(expr)
    return expr.collect(fcts)
    """
    return expr


def oprint(*args, dict_mode=False):
    """
    Debug printing for iterated (list/tuple/dict/set) objects. args is
    of form (title1, object1, title2, object2, ...) and prints:

        title1 = object1
        title2 = object2
        ...

    If you only wish to print a title set object = None.
    """

    if isinstance(args[0], str) or args[0] is None:
        titles = args[0::2]
        objs = args[1::2]
        strs = [
            ostr(obj, dict_mode) if obj is not None else None
            for obj in objs
        ]
        n = max((
            len(title)
            for title, s in zip(titles, strs)
            if s is not None and '\n' not in s
        ), default=0)

        for title, s in zip(titles, strs):
            if s is None:
                print(title)
            else:
                npad = n - len(title)
                if '\n' in s:
                    print(title + ':\n' + s)
                else:
                    print(title + npad * ' ' + ' = ' + s)
    else:
        for arg in args:
            print(ostr(arg, dict_mode))


class Eprint:

    ColorCode = {'black': '0;30', 'bright gray': '0;37',
                    'blue': '0;34', 'white': '1;37',
                    'green': '0;32', 'bright blue': '1;34',
                    'cyan': '0;36', 'bright green': '1;32',
                    'red': '0;31', 'bright cyan': '1;36',
                    'purple': '0;35', 'bright red': '1;31',
                    'yellow': '0;33', 'bright purple': '1;35',
                    'dark gray': '1;30', 'bright yellow': '1;33',
                    'normal': '0'}

    InvColorCode = dict(list(zip(list(ColorCode.values()), list(ColorCode.keys()))))

    normal = ''
    base = ''
    fct = ''
    deriv = ''
    bold = ''

    defaults = {('win', 'base'): 'blue', ('unix', 'base'): 'blue',
                ('win', 'fct'): 'red', ('unix', 'fct'): 'red',
                ('win', 'deriv'): 'cyan', ('unix', 'deriv'): 'cyan'}

    def __init__(self, base=None, fct=None, deriv=None, on=True, debug=False):
        if on:
            OS = 'unix'
            if 'win' in sys.platform and 'darwin' not in sys.platform:
                OS = 'win'

            if base is None:
                Eprint.base = Eprint.ColorCode[Eprint.defaults[(OS, 'base')]]
            else:
                Eprint.base = Eprint.ColorCode[base]
            if fct is None:
                Eprint.fct = Eprint.ColorCode[Eprint.defaults[(OS, 'fct')]]
            else:
                Eprint.fct = Eprint.ColorCode[fct]
            if deriv is None:
                Eprint.deriv = Eprint.ColorCode[Eprint.defaults[(OS, 'deriv')]]
            else:
                Eprint.deriv = Eprint.ColorCode[deriv]
            Eprint.normal = '\033[0m'

            if debug:
                print('Enhanced Printing is on:')
                print('Base/Blade color is ' + Eprint.InvColorCode[Eprint.base])
                print('Function color is ' + Eprint.InvColorCode[Eprint.fct])
                print('Derivative color is ' + Eprint.InvColorCode[Eprint.deriv] + '\n')

            Eprint.base = '\033[' + Eprint.base + 'm'
            Eprint.fct = '\033[' + Eprint.fct + 'm'
            Eprint.deriv = '\033[' + Eprint.deriv + 'm'

    @staticmethod
    def Base(s):
        return Eprint.base + s + Eprint.normal

    @staticmethod
    def Fct(s):
        return Eprint.fct + s + Eprint.normal

    @staticmethod
    def Deriv(s):
        return Eprint.deriv + s + Eprint.normal

    @staticmethod
    def Strip(s):
        new_s = s.replace(Eprint.base, '')
        new_s = new_s.replace(Eprint.normal, '')
        return new_s


class GaPrinter(StrPrinter):

    _default_settings = ChainMap({
        # if true, print dicts with `->` instead of `:`, one entry per line
        "dict_mode": False,
    }, StrPrinter._default_settings)

    function_names = ('acos', 'acosh', 'acot', 'acoth', 'arg', 'asin', 'asinh',
                      'atan', 'atan2', 'atanh', 'ceiling', 'conjugate', 'cos',
                      'cosh', 'cot', 'coth', 'exp', 'floor', 'im', 'log', 're',
                      'root', 'sin', 'sinh', 'sqrt', 'sign', 'tan', 'tanh', 'Abs')

    def _print_Function(self, expr):
        name = expr.func.__name__

        if expr.func.nargs is not None:
            if name in GaPrinter.function_names:
                return expr.func.__name__ + "(%s)" % self.stringify(expr.args, ", ")

        return Eprint.Fct("%s" % (name,))

    def _print_Derivative(self, expr):
        # Break the following to support both py 2 & 3
        # function, *diff_args = expr.args
        function = expr.args[0]
        diff_args = expr.args[1:]

        xi = []
        ni = []
        for x, n in diff_args:
            if x in xi:
                i = xi.index(x)
                ni[i] += n
            else:
                xi.append(self._print(x))
                ni.append(n)

        s = 'D'
        for x, n in zip(xi, ni):
            s += '{' + str(x) + '}'
            if n > 1:
                s += '^' + str(n)
        s += str(self._print(function))
        return Eprint.Deriv(s)

    def _print_dict(self, expr):
        if not self._settings['dict_mode']:
            return super()._print_dict(expr)

        return '\n'.join(
            '{} -> {}'.format(self._print(k), self._print(v))
            for k, v in expr.items()
        )


class GaPrintable:
    """ Mixin class providing default implementations of printing hooks """
    def __ga_print_str__(self):
        if GaLatexPrinter.latex_flg:
            return GaLatexPrinter().doprint(self)
        else:
            return GaPrinter().doprint(self)

    def __repr__(self):
        return GaPrinter().doprint(self)

    def _repr_latex_(self):
        # IPython expects latex in text mode, so we wrap in an environment
        latex_str = GaLatexPrinter().doprint(self)
        return r'\begin{equation*} ' + latex_str + r' \end{equation*}'


Basic.__ga_print_str__ = GaPrintable.__ga_print_str__
Matrix.__ga_print_str__ = GaPrintable.__ga_print_str__

Basic.__repr__ = GaPrintable.__repr__
Matrix.__repr__ = GaPrintable.__repr__


# This is the lesser of two evils. Previously, we overwrote `Basic.__str__` in
# order to customise `print(sympy)`. This broke a bunch of assumptions inside
# sympy, so isn't safe. Instead of clobbering `__str__`, we add a
# `__ga_print_str__` attribute, and have `print` use it if present.
_old_print = builtins.print


@functools.wraps(_old_print)
def _print(*values, **kwargs):
    values_new = []
    for v in values:
        try:
            f = type(v).__ga_print_str__
        except AttributeError:
            values_new.append(v)
        else:
            values_new.append(f(v))
    _old_print(*values_new, **kwargs)


builtins.print = _print


def enhance_print():
    Eprint()
    return


class GaLatexPrinter(LatexPrinter):
    r"""
    The latex printer is turned on with the function (in ga.py) -

        Format(Fmode=True, Dmode=True, ipy=False)

    where Fmode is the function printing mode that surpresses printing arguments,
    Dmode is the derivative printing mode that does not use fractions, and
    ipy=True is the Ipython notebook mode that does not redirect the print output.

    The latex output is post processed and displayed with the function (in GAPrint.py) -

        xpdf(filename='tmplatex.tex', debug=False)

    where filename is the name of the tex file one would keep for future
    inclusion in documents and debug=True would display the tex file
    immediately.

    There are three options for printing multivectors in latex.  They are
    acessed with the multivector member function -

        A.Fmt(self, fmt=1, title=None)

    where fmt=1, 2, or 3 determines whether the entire multivector A is
    printed entirely on one line, or one grade is printed per line, or
    one base is printed per line.  If title is not None then the latex
    string generated is of the form -

        title+' = '+str(A)

    where it is assumed that title is a latex math mode string. If title
    contains '%' it is treated as a pure latex math mode string.  If it
    does not contain '%' then the following character mappings are applied -

        'grad' -> '\bm{\nabla} '
        '*'    -> ''
        '^'    -> '\W '
        '|'    -> '\cdot '
        '>'    -> '\lfloor '
        '<'    -> '\rfloor '

    In the case of a print statement of the form -

        print title, A

    everthing in the title processing still applies except that the multivector
    formatting is one multivector per line.

    For print statements of the form -

        print title

    where no program variables are printed if title contains '#' then title
    is printed as regular latex line.  If title does not contain '#' then
    title is printed in equation mode. '%' has the same effect in title as
    in the Fmt() member function.
    """
    # overrides of base class settings, and new settings for our printers
    _default_settings = ChainMap({
        "mat_str": "array",
        "omit_function_args": False,
        "omit_partial_derivative_fraction": False,
    }, LatexPrinter._default_settings)

    latex_flg = False
    latex_str = ''
    ipy = False

    preamble = \
"""
\\pagestyle{empty}
\\usepackage[latin1]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{amsbsy}
\\usepackage{tensor}
\\usepackage{listings}
\\usepackage{color}
\\usepackage{xcolor}
\\usepackage{bm}
\\usepackage{breqn}
\\definecolor{gray}{rgb}{0.95,0.95,0.95}
\\setlength{\\parindent}{0pt}
\\DeclareMathOperator{\\Tr}{Tr}
\\DeclareMathOperator{\\Adj}{Adj}
\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}
\\newcommand{\\lp}{\\left (}
\\newcommand{\\rp}{\\right )}
\\newcommand{\\paren}[1]{\\lp {#1} \\rp}
\\newcommand{\\half}{\\frac{1}{2}}
\\newcommand{\\llt}{\\left <}
\\newcommand{\\rgt}{\\right >}
\\newcommand{\\abs}[1]{\\left |{#1}\\right | }
\\newcommand{\\pdiff}[2]{\\bfrac{\\partial {#1}}{\\partial {#2}}}
\\newcommand{\\lbrc}{\\left \\{}
\\newcommand{\\rbrc}{\\right \\}}
\\newcommand{\\W}{\\wedge}
\\newcommand{\\prm}[1]{{#1}'}
\\newcommand{\\ddt}[1]{\\bfrac{d{#1}}{dt}}
\\newcommand{\\R}{\\dagger}
\\newcommand{\\deriv}[3]{\\bfrac{d^{#3}#1}{d{#2}^{#3}}}
\\newcommand{\\grade}[1]{\\left < {#1} \\right >}
\\newcommand{\\f}[2]{{#1}\\lp{#2}\\rp}
\\newcommand{\\eval}[2]{\\left . {#1} \\right |_{#2}}
\\newcommand{\\Nabla}{\\boldsymbol{\\nabla}}
\\newcommand{\\eb}{\\boldsymbol{e}}
\\usepackage{float}
\\floatstyle{plain} % optionally change the style of the new float
\\newfloat{Code}{H}{myc}
\\lstloadlanguages{Python}

\\begin{document}
"""
    postscript = '\\end{document}\n'
    macros = '\\newcommand{\\f}[2]{{#1}\\left ({#2}\\right )}'

    # translate name, supers and subs to tex keywords
    greek = set(['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                 'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu',
                 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon',
                 'phi', 'chi', 'psi', 'omega'])

    greek_translated = {'lamda': 'lambda', 'Lamda': 'Lambda'}

    other = set(['aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth',
                 'hbar', 'hslash', 'mho', 'infty'])

    special_alphabet = list(reversed(sorted(list(greek) + list(other), key=len)))

    @staticmethod
    def redirect():
        GaLatexPrinter.latex_flg = True
        if GaLatexPrinter.ipy:
            pass
        else:
            GaLatexPrinter.stdout = sys.stdout
            sys.stdout = io.StringIO()

    @staticmethod
    def restore():
        if GaLatexPrinter.latex_flg:
            if not GaLatexPrinter.ipy:
                GaLatexPrinter.latex_str += sys.stdout.getvalue()
            GaLatexPrinter.latex_flg = False
            if not GaLatexPrinter.ipy:
                sys.stdout = GaLatexPrinter.stdout

    def _print_Pow(self, expr):
        base = self._print(expr.base)
        if ('_' in base or '^' in base) and 'cdot' not in base:
            mode = True
        else:
            mode = False

        # Treat x**Rational(1, n) as special case
        if expr.exp.is_Rational and abs(expr.exp.p) == 1 and expr.exp.q != 1:
            #base = self._print(expr.base)
            expq = expr.exp.q

            if expq == 2:
                tex = r"\sqrt{%s}" % base
            elif self._settings['itex']:
                tex = r"\root{%d}{%s}" % (expq, base)
            else:
                tex = r"\sqrt[%d]{%s}" % (expq, base)

            if expr.exp.is_negative:
                return r"\frac{1}{%s}" % tex
            else:
                return tex
        elif self._settings['fold_frac_powers'] \
            and expr.exp.is_Rational \
                and expr.exp.q != 1:
            base, p, q = self._print(expr.base), expr.exp.p, expr.exp.q
            if mode:
                return r"{\left ( %s \right )}^{%s/%s}" % (base, p, q)
            else:
                return r"%s^{%s/%s}" % (base, p, q)

        elif expr.exp.is_Rational and expr.exp.is_negative and expr.base.is_Function:
            # Things like 1/x
            return r"\frac{%s}{%s}" % \
                (1, self._print(Pow(expr.base, -expr.exp)))
        else:
            if expr.base.is_Function:
                return r"{%s}^{%s}" % (self._print(expr.base), self._print(expr.exp))
            else:
                if expr.is_commutative and expr.exp == -1:
                    #solves issue 1030
                    #As Mul always simplify 1/x to x**-1
                    #The objective is achieved with this hack
                    #first we get the latex for -1 * expr,
                    #which is a Mul expression
                    tex = self._print(S.NegativeOne * expr).strip()
                    #the result comes with a minus and a space, so we remove
                    if tex[:1] == "-":
                        return tex[1:].strip()
                if self._needs_brackets(expr.base):
                    tex = r"\left(%s\right)^{%s}"
                else:
                    if mode:
                        tex = r"{\left ( %s \right )}^{%s}"
                    else:
                        tex = r"%s^{%s}"

                return tex % (self._print(expr.base),
                              self._print(expr.exp))

    def _print_Symbol(self, expr, style='plain'):

        def str_symbol(name_str):

            def translate(s):
                tmp = s

                parse_dict = {}
                i_sub = 1

                for glyph in GaLatexPrinter.special_alphabet:
                    escaped_glyph = '\\' + glyph
                    if glyph in tmp:
                        parse_sym = '????' + str(i_sub)
                        i_sub += 1
                        # If this glyph is already escaped, avoid escaping again
                        translated_glyph = (escaped_glyph + ' ') if escaped_glyph not in tmp else glyph
                        parse_dict[parse_sym] = translated_glyph
                        tmp = tmp.replace(glyph, parse_sym)

                for parse_sym in parse_dict:
                    tmp = tmp.replace(parse_sym, parse_dict[parse_sym])

                for glyph in GaLatexPrinter.greek_translated:
                    if glyph in tmp:
                        tmp = tmp.replace(glyph, GaLatexPrinter.greek_translated[glyph])

                return tmp

            name, supers, subs = split_super_sub(name_str)

            name = translate(name)

            if style == 'bold':
                name = '\\boldsymbol{' + name + '}'

            supers = list(map(translate, supers))
            subs = list(map(translate, subs))

            # glue all items together:
            if len(supers) > 0:
                name += "^{%s}" % " ".join(supers)
            if len(subs) > 0:
                name += "_{%s}" % " ".join(subs)

            return name

        if expr in self._settings['symbol_names']:
            return self._settings['symbol_names'][expr]

        return str_symbol(expr.name)

    def _print_Function(self, expr, exp=None):

        func = expr.func.__name__
        name = func
        if hasattr(self, '_print_' + func):
            return getattr(self, '_print_' + func)(expr, exp)
        else:
            args = [str(self._print(arg)) for arg in expr.args]

            # How inverse trig functions should be displayed, formats are:
            # abbreviated: asin, full: arcsin, power: sin^-1
            inv_trig_style = self._settings['inv_trig_style']
            # If we are dealing with a power-style inverse trig function
            inv_trig_power_case = False
            # If it is applicable to fold the argument brackets
            can_fold_brackets = self._settings['fold_func_brackets'] and \
                len(args) == 1 and not self._needs_function_brackets(expr.args[0])

            inv_trig_table = ["asin", "acos", "atan", "acot", "acosh", "asinh", "atanh"]

            # If the function is an inverse trig function, handle the style
            if func in inv_trig_table:
                if inv_trig_style == "abbreviated":
                    func = func
                elif inv_trig_style == "full":
                    func = "arc" + func[1:]
                elif inv_trig_style == "power":
                    func = func[1:]
                    inv_trig_power_case = True

                    # Can never fold brackets if we're raised to a power
                    if exp is not None:
                        can_fold_brackets = False

            if inv_trig_power_case:
                if func in accepted_latex_functions:
                    name = r"\%s^{-1}" % func
                else:
                    name = r"\operatorname{%s}^{-1}" % func
            elif exp is not None:
                if func in accepted_latex_functions:
                    name = r"\%s^{%s}" % (func, exp)
                else:
                    name = latex(Symbol(func)) + ' '
                    if '_' in func or '^' in func:
                        name = r'{\left ( ' + name + r'\right ) }^{' + exp + '}'
                    else:
                        name += '^{' + exp + '}'
            else:
                if func in accepted_latex_functions:
                    name = r"\%s" % func
                else:
                    name = latex(Symbol(func)) + ' '
                    if exp is not None:
                        if '_' in name or '^' in name:
                            name = r'\left ( ' + name + r'\right )^{' + exp + '}'
                        else:
                            name += '^{' + exp + '}'

            if can_fold_brackets:
                if func in accepted_latex_functions:
                    # Wrap argument safely to avoid parse-time conflicts
                    # with the function name itself
                    name += r" {%s}"
                else:
                    if not self._settings["omit_function_args"]:
                        name += r"%s"
            else:
                if func in accepted_latex_functions or not self._settings["omit_function_args"]:
                    name += r"{\left (%s \right )}"

            if inv_trig_power_case and exp is not None:
                name += r"^{%s}" % exp

            if func in accepted_latex_functions or not self._settings["omit_function_args"]:
                if len(args) == 1:
                    name = name % args[0]
                else:
                    name = name % ",".join(args)

            return name

    def _print_Derivative(self, expr):
        dim = len(expr.variables)
        imax = 1
        if dim == 1:
            if self._settings["omit_partial_derivative_fraction"]:
                tex = r"\partial_{%s}" % self._print(expr.variables[0])
            else:
                tex = r"\frac{\partial}{\partial %s}" % self._print(expr.variables[0])
        else:
            multiplicity, i, tex = [], 1, ""
            current = expr.variables[0]
            for symbol in expr.variables[1:]:
                if symbol == current:
                    i = i + 1
                else:
                    multiplicity.append((current, i))
                    current, i = symbol, 1
            else:
                imax = max(imax, i)
                multiplicity.append((current, i))

            if self._settings["omit_partial_derivative_fraction"]:
                tex = ''
                for x, i in multiplicity:
                    if i == 1:
                        tex += r"\partial_{%s}" % (self._print(x),)
                    else:
                        tex += r"\partial^{%i}_{%s}" % (i, self._print(x))
            else:
                for x, i in multiplicity:
                    if i == 1:
                        tex += r"\partial %s" % self._print(x)
                    else:
                        tex += r"\partial^{%s} %s" % (i, self._print(x))
                tex = r"\frac{\partial^{%s}}{%s} " % (dim, tex)

        if isinstance(expr.expr, AssocOp):
            s = r"%s\left(%s\right)" % (tex, self._print(expr.expr))
        else:
            s = r"%s %s" % (tex, self._print(expr.expr))
        return s

    def _print_Determinant(self, expr):
        # sympy `uses |X|` by default, we want `det (X)`
        return r"\det\left ( {}\right )".format(self._print(expr.args[0]))

    @staticmethod
    def latex(expr, **settings):

        if not isinstance(expr, list):
            return GaLatexPrinter(settings).doprint(expr)
        else:
            s = '\\begin{align*}'
            for x in expr:
                s += '\n & ' + latex(x) + ' \\\\'
            s += '\n\\end{align*}'
            return s


def latex(expr, **settings):
    return GaLatexPrinter(settings).doprint(expr)


def print_latex(expr, **settings):
    """Prints LaTeX representation of the given expression."""
    print(latex(expr, **settings))


def Format(Fmode: bool = True, Dmode: bool = True, inverse='full'):
    r"""
    Turns on latex printing with configurable options.

    This redirects printer output so that latex compiler can capture it.

    ``Format()`` is also required for printing from *ipython notebook* (note that ``xpdf()`` is not needed to print from *ipython notebook*).

    Parameters
    ----------
    Fmode:
        * ``True`` -- Print functions without argument list, :math:`f`
        * ``False`` -- Print functions with standard *sympy* latex formatting, :math:`{{f}\lp {x,y,z} \rp }`
    Dmode:
        * ``True`` -- Print partial derivatives with condensed notation, :math:`\partial_{x}f`
        * ``False`` -- Print partial derivatives with standard *sympy* latex formatting, :math:`\pdiff{f}{x}`
    """
    global Format_cnt

    GaLatexPrinter.set_global_settings(
        omit_partial_derivative_fraction=Dmode,
        omit_function_args=Fmode,
        inv_trig_style=inverse,
    )

    if Format_cnt == 0:
        Format_cnt += 1

        GaLatexPrinter.latex_flg = True
        GaLatexPrinter.redirect()

        if isinteractive():
            init_printing(use_latex='mathjax')

    return


def _texify(s: str) -> str:
    """ Convert python GA operator notation to LaTeX """
    repl_pairs = [
        (r'\|', r'\cdot '),
        (r'\^(?!{)', r'\W '),
        (r'\*', ' '),
        (r'\brgrad\b', r'\bar{\boldsymbol{\nabla}} '),
        (r'\bgrad\b', r'\boldsymbol{\nabla} '),
        (r'>>', r' \times '),
        (r'<<', r' \bar{\times} '),
        (r'<', r'\rfloor '),
        (r'>', r'\lfloor '),
    ]

    def repl_func(m):
        # only one group will be present, use the corresponding match
        return next(
            r
            for (p, r), g in zip(repl_pairs, m.groups())
            if g is not None
        )
    pattern = '|'.join("({})".format(p) for p, _ in repl_pairs)
    return re.sub(pattern, repl_func, s)


def tex(paper=(14, 11), debug=False, prog=False, pt='10pt'):
    r"""
    Post processes LaTeX output (see comments below), adds preamble and
    postscript.

    This postprocessing has two main behaviors:

    1. Converting strings on the left hand side of the last ``=`` into TeX.
       This translates the ``*``, ``^``, ``|``, ``>``, ``<``, ``<<``, ``>>``,
       ``grad``, and ``rgrad`` operators of galgebra into the appropriate latex
       operators. If there is no ``=`` in the line, no conversion is applied.

    2. Wrapping lines of latex into ``equation*`` environments if they are not
       already in environments, and moving labels that were prepended outside
       ``align`` environments inside those environments.

    Both behaviors are applied line by line, unless a line starts with the
    following text:

    ``#%`` or ``%``
        Disables only behavior 1 for the rest of the line.

    ``##``
        Disables behaviors 1 and 2 until the end of the next line starting with
        ``##``. This includes processing any of the other special characters,
        which will be emitted verbatim.

    ``#``
        Disables behaviors 1 and 2 for the rest of the line.

    We assume that if :func:`tex` is called, then :func:`Format` has been called
    at the beginning of the program.
    """

    latex_str = GaLatexPrinter.latex_str + sys.stdout.getvalue()
    GaLatexPrinter.latex_str = ''
    GaLatexPrinter.restore()
    r"""
    Each line in the latex_str is interpreted to be an equation or align
    environment.  If the line does not begin with '\begin{align*}' then
    'begin{equation*}' will be added to the beginning of the line and
    '\end{equation*}' to the end of the line.
    The latex strings generated by galgebra and sympy expressions for
    printing must not contain '\n' except as the final character.  Thus
    all '\n' must be removed from a compound (not a simple type) expression
    and a '\n' added to the end of the string to delimit it when the string
    is generated.
    """
    latex_lst = latex_str.split('\n')
    latex_str = ''

    code_flg = False

    for latex_line in latex_lst:
        if not latex_line:
            pass
        elif latex_line.startswith('##'):
            # a post-processing toggle used by `Print_Function`
            code_flg = not code_flg
            latex_line = latex_line[2:]
        elif code_flg:
            pass
        elif latex_line.startswith('#') and not latex_line.startswith('#%'):
            # do not process this line
            latex_line = latex_line[1:]
        else:
            # two different spellings of "do not process the LHS"
            if latex_line.startswith('%'):
                latex_line = latex_line[1:]
            elif latex_line.startswith('#%'):
                latex_line = latex_line[2:]
            # otherwise, process it if we can find it
            elif '=' in latex_line:
                lhs, latex_line = latex_line.rsplit('=', 1)
                latex_line = _texify(lhs) + '=' + latex_line

            # in either case, perform the environment wrapping
            if r'\begin{align*}' in latex_line:
                latex_line = r'\begin{align*} ' + latex_line.replace(r'\begin{align*}', '', 1).lstrip()
            else:
                latex_line = r'\begin{equation*} ' + latex_line.strip() + r' \end{equation*}'

        latex_str += latex_line + '\n'

    latex_str = latex_str.replace('\n\n', '\n')

    if prog:
        with open(sys.argv[0], 'r') as prog_file:
            prog_str = prog_file.read()
        prog_str = '{\\Large \\bf Program:}\\begin{lstlisting}[language=Python,showspaces=false,' + \
                   'showstringspaces=false]\n' + \
                   prog_str + '\n\\end{lstlisting}\n {\\Large \\bf Code Output:} \n'
        latex_str = prog_str + latex_str

    if debug:
        print(latex_str)

    if paper == 'letter':
        paper_size = \
"""
\\documentclass[@10pt@,fleqn]{report}
"""
    else:
        paper_size = \
"""
\\documentclass[@10pt@,fleqn]{report}
\\usepackage[vcentering]{geometry}
"""
        if paper == 'landscape':
            paper = [11, 8.5]
        paper_size += '\\geometry{papersize={' + str(paper[0]) + \
                      'in,' + str(paper[1]) + 'in},total={' + str(paper[0] - 1) + \
                      'in,' + str(paper[1] - 1) + 'in}}\n'

    paper_size = paper_size.replace('@10pt@', pt)
    latex_str = paper_size + GaLatexPrinter.preamble + latex_str + GaLatexPrinter.postscript

    return latex_str


def xpdf(filename=None, paper=(14, 11), crop=False, png=False, prog=False, debug=False, pt='10pt', pdfprog='pdflatex'):

    """
    Post processes LaTeX output (see comments below), adds preamble and
    postscript, generates tex file, inputs file to latex, displays resulting
    pdf file.

    Arg         Value       Result
    pdfprog    'pdflatex'   Use pdfprog to generate pdf output, only generate tex if pdfprog is None
    crop        True        Use "pdfcrop" to crop output file (pdfcrop must be installed, linux only)
    png         True        Use "convert" to produce png output (imagemagick must be installed, linux only)

    We assume that if xpdf() is called then Format() has been called at the beginning of the program.
    """

    sys_cmd = SYS_CMD[sys.platform]

    latex_str = tex(paper=paper, debug=debug, prog=prog, pt=pt)

    if filename is None:
        pyfilename = sys.argv[0]
        rootfilename = pyfilename.replace('.py', '')
        filename = rootfilename + '.tex'

    if debug:
        print('latex file =', filename)

    latex_file = open(filename, 'w')
    latex_file.write(latex_str)
    latex_file.close()

    latex_str = None

    if pdfprog is None:
        return

    pdflatex = shutil.which(pdfprog)

    if debug:
        print('pdflatex path =', pdflatex)

    if pdfprog is not None:
        if debug:  # Display latex excution output for debugging purposes
            os.system(pdfprog + ' ' + filename[:-4])
        else:  # Works for Linux don't know about Windows
            os.system(pdfprog + ' ' + filename[:-4] + sys_cmd['null'])

        print_cmd = sys_cmd['evince'] + ' ' + filename[:-4] + '.pdf ' + sys_cmd['&']
        print(print_cmd)

        os.system(print_cmd)
        eval(input('!!!!Return to continue!!!!\n'))

        if debug:
            os.system(sys_cmd['rm'] + ' ' + filename[:-4] + '.aux ' + filename[:-4] + '.log')
        else:
            os.system(sys_cmd['rm'] + ' ' + filename[:-4] + '.aux ' + filename[:-4] + '.log ' + filename[:-4] + '.tex')
        if crop:
            os.system('pdfcrop ' + filename[:-4] + '.pdf')
            os.remove(filename[:-4] + '.pdf')
            os.rename(filename[:-4] + '-crop.pdf', filename[:-4] + '.pdf')
        if png:
            os.system('Pdf2Png ' + filename[:-4])
    return


def xdvi(filename=None, debug=False, paper=(14, 11)):
    xpdf(filename=filename, paper=paper, crop=False, png=False, prog=False, debug=debug, pt='10pt')
    return


def LatexFormat(Fmode=True, Dmode=True, ipy=False):
    GaLatexPrinter.set_global_settings(
        omit_partial_derivative_fraction=Dmode,
        omit_function_args=Fmode
    )
    GaLatexPrinter.ipy = ipy
    GaLatexPrinter.redirect()
    return


off_mode = False


def Get_Program(off=False):
    global off_mode
    off_mode = off


def Print_Function():
    if off_mode:
        return

    tmp_str = inspect.getsource(inspect.currentframe().f_back)
    if GaLatexPrinter.latex_flg:
        #print '#Code for '+fct_name
        print(r'##\begin{lstlisting}[language=Python,showspaces=false,'
              r'showstringspaces=false,backgroundcolor=\color{gray},frame=single]')
        print(tmp_str)
        print('##\\end{lstlisting}')
        print('#Code Output:')
    else:
        print('\n' + 80 * '*')
        #print '\nCode for '+fct_name
        print(tmp_str)
        print('Code output:\n')
    return


_eval_global_dict = {}
_eval_parse_order = []


def def_prec(gd: dict, op_ord: str = '<>|,^,*') -> None:
    """
    This is used with the ``GAeval()`` function to evaluate a string representing a multivector expression with a revised operator precedence.

    Parameters
    ----------
    gd :
        The ``globals()`` dictionary to lookup variable names in.
    op_ord :
        The order of operator precedence from high to low with groups of equal precedence separated by commas.
        The default precedence, ``'<>|,^,*'``, is that used by Hestenes (:cite:`Hestenes`, p7, :cite:`Doran`, p38).
        This means that the ``<``, ``>``, and ``|`` operations have equal
        precedence, followed by ``^``, and lastly ``*``.
    """
    global _eval_global_dict, _eval_parse_order
    op_ord_list = op_ord.split(',')
    _parser.validate_op_order(op_ord_list)
    _eval_global_dict = gd
    _eval_parse_order = op_ord_list


def GAeval(s: str, pstr: bool = False):
    """
    Evaluate a multivector expression string ``s``.

    The operator precedence and variable values within the string are
    controlled by :func:`def_prec`. The documentation for that function
    describes the default precedence.

    The implementation works by adding parenthesis to the input string ``s``
    according to the requested precedence, and then calling :func:`eval` on the
    result.

    For example consider where ``X``, ``Y``, ``Z``, and ``W`` are multivectors::

        def_prec(globals())
        V = GAeval('X|Y^Z*W')

    The *sympy* variable ``V`` would evaluate to ``((X|Y)^Z)*W``.

    Parameters
    ----------
    s :
        The string to evaluate.
    pstr :
        If ``True``, the values of ``s`` and ``s`` with parenthesis added to
        enforce operator precedence are printed.
    """

    seval = _parser.parse_line(s, _eval_parse_order)
    if pstr:
        print(s)
        print(seval)
    return eval(seval, _eval_global_dict)


def Fmt(obj, fmt=0):
    if isinstance(obj, (list, tuple, dict)):
        n = len(obj)
        if isinstance(obj, list):
            ldelim = '['
            rdelim = ']'
        elif isinstance(obj, dict):
            ldelim = r'\{'
            rdelim = r'\}'
        else:
            ldelim = '('
            rdelim = ')'
        if fmt == 1:
            latex_str = r' \left ' + ldelim + r' \begin{array}{' + n*'c' + '} '
            for cell in obj:
                if isinstance(obj, dict):
                    #cell.title = None
                    latex_cell = latex(cell) + ' : ' + latex(obj[cell])
                else:
                    #title = cell.title
                    #cell.title = None
                    latex_cell = latex(cell)
                latex_cell = latex_cell.replace('\n', ' ')
                latex_str += latex_cell + ', & '
                #cell.title = title
            latex_str = latex_str[:-4]
            latex_str += r'\\ \end{array} \right ' + rdelim + ' \n'
        else:
            latex_str = ''
            i = 1
            for cell in obj:
                #title = cell.title
                #cell.title = None
                latex_cell = latex(cell)
                latex_cell = latex_cell.replace('\n', ' ')
                #cell.title = title
                if i == 1:
                    latex_str += r'\begin{array}{c} \left ' + ldelim + r' ' + latex_cell + r', \right. \\ '
                elif i == n:
                    latex_str += r' \left. ' + latex_cell + r'\right ' + rdelim + r' \\ \end{array}'
                else:
                    latex_str += r' ' + latex_cell + r', \\'
                i += 1
        if isinteractive():  # For Ipython notebook
            latex_str = r'\begin{equation*} ' + latex_str + r'\end{equation*}'
            return latex_str
        else:
            return latex_str

    elif isinstance(obj, int):
        LatexPrinter.set_global_settings(galgebra_mv_fmt=obj)
        return
    else:
        raise TypeError(str(type(obj)) + ' not allowed arg type in Fmt')


class _WithSettings(GaPrintable):
    """ Helper class to attach print settings to an object """
    def __init__(self, obj, settings: dict = {}):
        self._obj = obj
        self._settings = settings

    def __do_print(self, printer):
        # make a copy of the printer with the specified setting applied
        new_printer = copy.copy(printer)
        new_printer._settings = copy.copy(new_printer._settings)
        new_printer._settings.update(self._settings)
        return new_printer.doprint(self._obj)

    _latex = _pretty = _sympystr = __do_print


class _FmtResult(GaPrintable):
    """ Object returned from .Fmt methods, which can be printed as latex """
    def __new__(cls, obj, label: str) -> GaPrintable:
        if label is None:
            return obj
        self = super().__new__(cls)
        self._obj = obj
        self._label = label
        return self

    def _latex(self, printer):
        return self._label + ' = ' + printer.doprint(self._obj)

    def _sympystr(self, printer):
        return self._label + ' = ' + printer.doprint(self._obj)
