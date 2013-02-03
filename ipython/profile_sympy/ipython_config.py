import sys,os

c = get_config()
app = c.InteractiveShellApp
c.InteractiveShell.editor = 'geany'

# This can be used at any point in a config file to load a sub config
# and merge it into the current one.
load_subconfig('ipython_config.py', profile='default')

lines = """from sympy import *
from sympy.galgebra.GA import *
import sympy.galgebra.latex_ex as tex"""

# You have to make sure that attributes that are containers already
# exist before using them.  Simple assigning a new list will override
# all previous values.

if hasattr(app, 'exec_lines'):
    app.exec_lines.append(lines)
else:
    app.exec_lines = [lines]

# Load the sympy_printing extension to enable nice printing of sympy expr's.
if hasattr(app, 'extensions'):
    app.extensions.append('GAprinting')
else:
    app.extensions = ['GAprinting']


if hasattr(app, 'extensions'):
    app.extensions.append('laga_ops')
else:
    app.extensions = ['laga_ops']

