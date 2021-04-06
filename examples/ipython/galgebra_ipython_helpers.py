import subprocess
import sys
import os
from IPython.display import display_pretty


def check(name):
    run(name)

    with open(name + '.tex', 'r') as f:
        # can't use display.Latex here, it would result in CSS comparisons in the output.
        # using `display` forces this to be a separate output to any stdout from above.
        display_pretty(f.read(), raw=True)


def run(name):
    # this makes Python < 3.9 behave like 3.9
    abs_name = os.path.abspath(name)
    # stdout and stderr do not seem to go the right place in jupyter
    p = subprocess.run(
        [sys.executable, '-Wdefault', abs_name + '.py'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    sys.stdout.write(p.stdout)
    sys.stdout.flush()
    # remove the absolute paths from deprecation warnings
    sys.stderr.write(p.stderr.replace(abs_name, name))
    sys.stderr.flush()

    # this makes the error easier to read in nbval
    if p.returncode:
        raise RuntimeError("The script raised an exception:\n\n" + p.stderr)
