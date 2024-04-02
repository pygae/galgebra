# Author:  Greg Grunberg
# Last revised:  2023-01-13

# INSTRUCTIONS:
#    
# Place this GAlgebrainit.py file in the "galgebra" folder, which should
# already contain the other GAlgebra modules 
#
#     atoms.py          dop.py            ga.py
#     lt.py             metric.py         mv.py
#     printer.py        gprinter.py       utils.py .
#
# To invoke this program within a Jupyter notebok, write the command 
#
#     from galgebra.GAlgebraInit import *
#
# in the notebook's first In[ ] cell and execute that cell.

# Make SymPy available to this program
import sympy 
from sympy import *

# Make GAlgebra available to this program.
import galgebra
from galgebra.ga import *  
from galgebra.mv import *
from galgebra.lt import *
from galgebra.dop import *
from galgebra.printer import Fmt, GaPrinter, Format
    # Fmt:  sets display mode of a multivector's basis expansion.
    # GaPrinter:  makes GA output a little more readable.
    # Format:  turns on latex printer.
from galgebra.gprinter import gFormat, gprint
gFormat()
    # Default `Fmode=True` suppresses display of the arguments of 
    #   multivector fields.  
    # Default `Dmode=True` causes partial differentiation
    #   operators to be displayed in shortened form.
Ga.dual_mode('Iinv+')  
    # Sets multivector dualization to be right multiplication by the 
    # by the inverse unit pseudoscalar (the convention used in the 
    # textbooks LAGA, VAGC, and GACS).

initializations_list = \
r"""\textsf{The following initialization commands were executed:}\\
\quad\texttt{from sys import version}\\
\quad\texttt{import sympy}\\
\quad\texttt{from sympy import *}\\
\quad\texttt{import galgebra}\\
\quad\texttt{from galgebra.ga import *}\\  
\quad\texttt{from galgebra.mv import *}\\
\quad\texttt{from galgebra.lt import *}\\   
\quad\texttt{from galgebra.dop import *}\\
\quad\texttt{from galgebra.printer import Fmt, GaPrinter, Format}\\
\quad\texttt{from galgebra.gprinter import gFormat, gprint}\\
\quad\texttt{gFormat()}\\~~~~~~~~
\quad\textsf{# default 'Fmode=True' suppresses arguments of multivector fields}\\~~~~~~~~
\quad\textsf{# default 'Dmode=True' displays partial differentiation operators in short form}\\
\quad\texttt{Ga.dual_mode('Iinv+')}\\~~~~~~~~
\quad\textsf{# dual and undual defined by }\mathbf{M}^\star = \mathbf{M}\mathbf{I}^{-1} \textsf{ and } \mathbf{M}^{-\star} = \mathbf{M}\mathbf{I}\\"""
gprint(initializations_list)

from sys import version
# Display versions of softwared being used.
gprint(r'\textsf{This notebook is now using} \\',
       r'\qquad\bullet~ \textsf{Python }', version[:5],
       r'\qquad\bullet~ \textsf{SymPy }', sympy.__version__[:8],
       r'\qquad\bullet~ \textsf{GAlgebra }', galgebra.__version__[:], r'.')
