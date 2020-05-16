import subprocess
import sys
from IPython.display import display_pretty


def check(name):
    run(name)

    with open(name + '.tex', 'r') as f:
        # can't use display.Latex here, it would result in CSS comparisons in the output.
        # using `display` forces this to be a separate output to any stdout from above.
        display_pretty(f.read(), raw=True)


def run(name):
    # stdout and stderr do not seem to go the right place in jupyter
    p = subprocess.run(
        [sys.executable, '-Wdefault', name + '.py'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    sys.stdout.write(p.stdout)
    sys.stdout.flush()
    sys.stderr.write(p.stderr)
    sys.stderr.flush()

    # this makes the error easier to read in nbval
    if p.returncode:
        raise RuntimeError("The script raised an exception:\n\n" + p.stderr)
