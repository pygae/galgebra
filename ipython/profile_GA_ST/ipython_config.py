import sys

GA_ST = """from sympy import *
GA import *
(t,x,y,z) = symbols('t x y z')
(et,ex,ey,ez) = MV.setup('e_t e_x e_y e_z',metric='[1,-1,-1,-1]',coords=(t,x,y,z))
Format(ipy=True)
print 'Spacetime mode'
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
    app.exec_lines.append(GA_ST)
else:
    app.exec_lines = [GA_ST]









