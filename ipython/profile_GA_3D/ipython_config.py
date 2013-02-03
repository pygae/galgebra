import sys

GA_3D = """from sympy import *
from GA import *
(x,y,z) = symbols('x y z')
(ex,ey,ez) = MV.setup('e_x e_y e_z',metric='[1,1,1]',coords=(x,y,z))
Format(ipy=True)
print '3D mode'
%load_ext GAprinting
"""

c = get_config()
app = c.InteractiveShellApp
c.InteractiveShell.editor = 'geany'

# This can be used at any point in a config file to load a sub config
# and merge it into the current one.
load_subconfig('ipython_config.py', profile='default')

# You have to make sure that attributes that are containers already
# exist before using them.  Simple assigning a new list will override
# all previous values.

if hasattr(app, 'exec_lines'):
    app.exec_lines.append(GA_3D)
else:
    app.exec_lines = [GA_3D]









