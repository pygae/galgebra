"""
A print function that pretty prints sympy Basic objects.

:moduleauthor: Brian Granger

Usage
=====

Once the extension is loaded, Sympy Basic objects are automatically
pretty-printed.

"""
#-----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from IPython.lib.latextools import latex_to_png
from IPython.testing import decorators as dec
# use @dec.skipif_not_sympy to skip tests requiring sympy

try:
    #from sympy import pretty, latex
    from sympy import pretty
    from GAPrint import latex
    from GA import MV
    from laga import parse_line
except ImportError:
    pass


#-----------------------------------------------------------------------------
# Definitions of magic functions for use with IPython
#-----------------------------------------------------------------------------

def print_basic_unicode(o, p, cycle):
    """A function to pretty print sympy Basic objects."""

    if cycle:
        return p.text('Basic(...)')
    if isinstance(o,MV):
        out = MV.str_rep(o)
    else:
        out = pretty(o, use_unicode=True)
    if '\n' in out:
        p.text(u'\n')
    p.text(out)


def print_png(o):
    """A function to display sympy expression using LaTex -> PNG."""
    #s = latex(o, mode='inline')
    # mathtext does not understand certain latex flags, so we try to replace
    # them with suitable subs.
    s = LaTeX(o,itex=True)
    s = s.replace('\\operatorname','')
    s = s.replace('\\overline', '\\bar')
    png = latex_to_png(s)
    return png

newcommands =  '$\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}$\n'+\
               '$\\newcommand{\\lp}{\\left (}$\n'+\
               '$\\newcommand{\\rp}{\\right )}$\n'+\
               '$\\newcommand{\\half}{\\frac{1}{2}}$\n'+\
               '$\\newcommand{\\llt}{\\left <}$\n'+\
               '$\\newcommand{\\rgt}{\\right >}$\n'+\
               '$\\newcommand{\\abs}[1]{\\left |{#1}\\right | }$\n'+\
               '$\\newcommand{\\pdiff}[2]{\\bfrac{\\partial {#1}}{\\partial {#2}}}$\n'+\
               '$\\newcommand{\\lbrc}{\\left \\{}$\n'+\
               '$\\newcommand{\\rbrc}{\\right \\}}$\n'+\
               '$\\newcommand{\\W}{\\wedge}$\n'+\
               "$\\newcommand{\\prm}[1]{{#1}'}$\n"+\
               '$\\newcommand{\\ddt}[1]{\\bfrac{d{#1}}{dt}}$\n'+\
               '$\\newcommand{\\R}{\\dagger}$\n'+\
               '$\\newcommand{\\bm}{\\boldsymbol}$\n'

_first_pass = True

def print_latex(o):
    global _first_pass,newcommands
    """A function to generate the latex representation of sympy expressions."""
    #s = latex(o, mode='equation', itex=True)
    s = latex(o,inline=False,itex=True)
    #s = s.replace('\\dag','\\dagger')
    if _first_pass:
        _first_pass = False
        s = newcommands+s
    return s

def parse_operators(line):
    return(line)

class QTransformer(object):
    # XXX: inheriting from PrefilterTransformer as documented gives TypeErrors,
    # but apparently is not needed after all
    priority = 99
    enabled = True
    def transform(self, line, continue_prompt):
        if line[0] == '>':
            line = parse_line(line[1:])
        return line

q_transformer = QTransformer()

_loaded = False

def load_ipython_extension(ip):
    """Load the extension in IPython."""
    global _loaded,q_transformer
    if not _loaded:
        plaintext_formatter = ip.display_formatter.formatters['text/plain']
        ip.prefilter_manager.register_transformer(q_transformer)

        for cls in (object, tuple, list, set, frozenset, dict, str):
            plaintext_formatter.for_type(cls, print_basic_unicode)

        plaintext_formatter.for_type_by_name(
            'sympy.core.basic', 'Basic', print_basic_unicode
        )
        plaintext_formatter.for_type_by_name(
            'sympy.matrices.matrices', 'Matrix', print_basic_unicode
        )
        plaintext_formatter.for_type_by_name(
            'GA', 'MV', print_basic_unicode
        )
        png_formatter = ip.display_formatter.formatters['image/png']

        png_formatter.for_type_by_name(
            'sympy.core.basic', 'Basic', print_png
        )

        latex_formatter = ip.display_formatter.formatters['text/latex']
        latex_formatter.for_type_by_name(
            'sympy.core.basic', 'Basic', print_latex
        )
        latex_formatter.for_type_by_name(
            'sympy.matrices.matrices', 'Matrix', print_latex
        )
        latex_formatter.for_type_by_name(
            'GA', 'MV', print_latex
        )
        _loaded = True
        print 'GA options loaded'

def unload_ipython_extension(ip):
    global q_transformer
    ip.prefilter_manager.unregister_transformer(q_transformer)
