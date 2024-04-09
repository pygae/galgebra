from .printer import GaLatexPrinter, isinteractive, Format_cnt, latex
import subprocess
import sys
import shutil

from sympy import init_printing

try:
    from IPython.display import display, Math
except ImportError:
    pass

SYS_CMD = {'linux2': {'rm': 'rm', 'evince': 'evince', 'null': ' > /dev/null', '&': '&'},
           'linux': {'rm': 'rm', 'evince': 'evince', 'null': ' > /dev/null', '&': '&'},
           'win32': {'rm': 'del', 'evince': 'start', 'null': ' > NUL', '&': ''},
           'darwin': {'rm': 'rm', 'evince': 'open', 'null': ' > /dev/null', '&': '&'}}


class LaTeX:
    # LaTeX data
    line_sep = """
************************************************************************
"""
    latex_flg = False
    latex_str = ''

    latex_preamble = """
\\pagestyle{empty}
\\usepackage[utf8]{inputenc}
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
\\newcommand{\\prm}[1]{{#1}^{\\prime}}
\\newcommand{\\ddt}[1]{\\bfrac{d{#1}}{dt}}
\\newcommand{\\R}{\\dagger}
\\newcommand{\\deriv}[3]{\\bfrac{d^{#3}#1}{d{#2}^{#3}}}
\\newcommand{\\grade}[2]{\\left < {#1} \\right >_{#2}}
\\newcommand{\\f}[2]{{#1}\\lp{#2}\\rp}
\\newcommand{\\eval}[2]{\\left . {#1} \\right |_{#2}}
\\newcommand{\\Nabla}{\\boldsymbol{\\nabla}}
\\newcommand{\\eb}{\\boldsymbol{e}}
\\newcommand{\\bs}[1]{\\boldsymbol{#1}}
\\newcommand{\\grad}{\\bs{\\nabla}}
\\usepackage{float}
\\floatstyle{plain} % optionally change the style of the new float
\\newfloat{Code}{H}{myc}
\\lstloadlanguages{Python}
\\begin{document}

"""

    ip_cmds = \
        [r'$$\DeclareMathOperator{\Tr}{Tr}$$',
            r'$$\DeclareMathOperator{\Adj}{Adj}$$',
            r'$$\newcommand{\bfrac}[2]{\displaystyle\frac{#1}{#2}}$$',
            r'$$\newcommand{\lp}{\left (}$$',
            r'$$\newcommand{\rp}{\right )}$$',
            r'$$\newcommand{\paren}[1]{\lp {#1} \rp}$$',
            r'$$\newcommand{\half}{\frac{1}{2}}$$',
            r'$$\newcommand{\llt}{\left <}$$',
            r'$$\newcommand{\rgt}{\right >}$$',
            r'$$\newcommand{\abs}[1]{\left |{#1}\right | }$$',
            r'$$\newcommand{\pdiff}[2]{\bfrac{\partial {#1}}{\partial {#2}}}$$',
            r'$$\newcommand{\npdiff}[3]{\bfrac{\partial^{#3} {#1}}{\partial {#2}^{#3}}}$$',
            r'$$\newcommand{\lbrc}{\left \{}$$',
            r'$$\newcommand{\rbrc}{\right \}}$$',
            r'$$\newcommand{\W}{\wedge}$$',
            r'$$\newcommand{\prm}[1]{{#1}^{\prime}}$$',
            r'$$\newcommand{\ddt}[1]{\bfrac{d{#1}}{dt}}$$',
            r'$$\newcommand{\R}{\dagger}$$',
            r'$$\newcommand{\deriv}[3]{\bfrac{d^{#3}#1}{d{#2}^{#3}}}$$',
            r'$$\newcommand{\grade}[2]{\left < {#1} \right >_{#2}}$$',
            r'$$\newcommand{\f}[2]{{#1}\lp {#2} \rp}$$',
            r'$$\newcommand{\eval}[2]{\left . {#1} \right |_{#2}}$$',
            r'$$\newcommand{\bs}[1]{\boldsymbol{#1}}$$',
            r'$$\newcommand{\grad}{\bs{\nabla}}$$']

# ***********************************************************************


def gFormat(Fmode: bool = True, Dmode: bool = True, inverse='full'):
    r"""
    Turns on latex printing with configurable options.

    This redirects printer output so that latex compiler can capture it.

    ``Format()`` is also required for printing from *ipython notebook* (note that ``xpdf()`` is not needed to print from *ipython notebook*).

    Parameters
    ----------
    Fmode:
        Value for the ``omit_function_args`` setting of
        :class:`GaLatexPrinter`.
    Dmode:
        Value for the ``omit_partial_derivative_fraction`` setting of
        :class:`GaLatexPrinter`.
    """
    global Format_cnt

    GaLatexPrinter.set_global_settings(
        omit_partial_derivative_fraction=Dmode,
        omit_function_args=Fmode,
        inv_trig_style=inverse,
    )

    if Format_cnt == 0:
        Format_cnt += 1

        LaTeX.latex_flg = True

        if isinteractive():
            init_printing(use_latex='mathjax')
            from IPython.display import Math, display
            cmds = '\n'.join(LaTeX.ip_cmds)
            display(Math(cmds))

    return


def gprint(*xargs):
    """
    Print latex or text from python script or latex from Jupyter Notebook/Lab

    """
    x = []
    fstr = ''
    new_eq_flg = False
    i = 0
    for xi in xargs:
        if isinstance(xi, str):
            if r'\\' in xi and i > 0:
                if isinteractive():
                    xi_rep = xi.replace(r'\\', r'\end{equation*}@\begin{equation*} ')
                else:
                    xi_rep = xi.replace(r'\\', r'\end{equation*}'+'\n'+r'\begin{equation*} ')
                new_eq_flg = True
                fstr += xi_rep
            else:
                fstr += xi
        elif isinstance(xi, type):
            if LaTeX.latex_flg:
                fstr += r' \text{'+str(xi)+'} '
            else:
                fstr += str(xi)
        else:
            if LaTeX.latex_flg:
                x.append(latex(xi))
                if new_eq_flg:
                    new_eq_flg = False
                fstr += r' %s '
            else:
                x.append(str(xi))
                fstr += r' %s '

        i += 1

    if LaTeX.latex_flg:
        if isinteractive():
            lstr = fstr % tuple(x)
            if '@' in lstr:
                lines = lstr.split('@')
                lines[0] = r'\begin{equation*} '+lines[0]
                lines[-1] += r'\end{equation*}'
                for line in lines:
                    display(Math(line))
            else:
                display(Math(lstr))
        else:
            LaTeX.latex_str += r'\begin{equation*} ' + (fstr % tuple(x)) + r'\end{equation*} '+'\n'
    else:
        print(fstr % tuple(x))

    return


def gxpdf(filename=None, paper=(14, 11), crop=False, png=False, prog=False, debug=False, pt='10pt', pdfprog='pdflatex', evince=True, rm=True, null=True):

    """
    Post processes LaTeX output (see comments below), adds preamble and
    postscript, generates tex file, inputs file to latex, displays resulting
    pdf file.

    Arg         Value       Result
    pdfprog    'pdflatex'   Use pdfprog to generate pdf output, only generate tex if pdfprog is None
    crop        True        Use "pdfcrop" to crop output file (pdfcrop must be installed, linux only)
    png         True        Use "convert" to produce png output (imagemagick must be installed, linux only)

    We assume that if gxpdf() is called then gFormat() has been called at the beginning of the program.
    """

    latex_str = paper_format(paper, pt)+LaTeX.latex_preamble+LaTeX.latex_str+r'\end{document}'

    if filename is None:
        pyfilename = sys.argv[0]
        rootfilename = pyfilename.replace('.py', '')
        tex_filename = rootfilename + '.tex'
        pdf_filename = rootfilename + '.pdf'
    else:
        tex_filename = filename
        pdf_filename = tex_filename.replace('.tex', '.pdf')
        rootfilename = tex_filename.replace('.tex', '')

    if debug:
        print('latex file =', filename)

    with open(tex_filename, 'w') as latex_file:
        latex_file.write(latex_str)

    if pdfprog is not None:
        sys_cmd = SYS_CMD[sys.platform]
        pdflatex = shutil.which(pdfprog)

        if debug:  # Display latex excution output for debugging purposes
            print('pdflatex path =', pdflatex)
            # os.system(pdfprog + ' ' + filename[:-4])
        else:  # Works for Linux don't know about Windows
            if null:
                subprocess.call([pdfprog, tex_filename, sys_cmd['null']])
            else:
                subprocess.call([pdfprog, tex_filename])
            # os.system(pdfprog + ' ' + filename[:-4] + sys_cmd['null'])

        if evince:
            subprocess.call([sys_cmd['evince'], pdf_filename])

        # eval(input('!!!!Return to continue!!!!\n'))

        if rm:
            if debug:
                subprocess.call([sys_cmd['rm'], rootfilename+'.aux ', rootfilename+'.log'])
            else:
                subprocess.call([sys_cmd['rm'], rootfilename+'.aux ', rootfilename+'.log ', rootfilename+'.tex'])
        if crop:
            subprocess.call(['pdfcrop', pdf_filename])
            subprocess.call(['rm', pdf_filename])
            subprocess.call(['mv', rootfilename+'-crop.pdf', pdf_filename])
        if png:
            subprocess.call(['Pdf2Png', rootfilename])
    return


def paper_format(paper, pt):  # Set size of paper and font size

    if paper == 'letter':
        paper_size = """
\\documentclass[@10pt@,fleqn]{book}
"""
    else:
        paper_size = """
\\documentclass[@10pt@,fleqn]{book}
\\usepackage[vcentering]{geometry}
"""
        if paper == 'landscape':
            paper = [11, 8.5]
        paper_size += '\\geometry{papersize={' + str(paper[0]) + \
                      'in,' + str(paper[1]) + 'in},total={' + str(paper[0] - 1) + \
                      'in,' + str(paper[1] - 1) + 'in}}\n'

    paper_size = paper_size.replace('@10pt@', pt)

    return paper_size
