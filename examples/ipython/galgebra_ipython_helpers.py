import subprocess
import sys
from IPython.display import display_pretty


def check(name):
    # todo: use subprocess.run once we drop Python 2.7
    p = subprocess.Popen([sys.executable, name + '.py'], stderr=subprocess.PIPE, universal_newlines=True)
    try:
        stdout, stderr = p.communicate()
    except:
        p.kill()
    if p.poll():
        raise RuntimeError("The script raised an exception:\n\n" + stderr)

    with open(name + '.tex', 'r') as f:
        # can't use display.Latex here, it would result in CSS comparisons in the output.

        # using `display` forces this to be a separate output to any stdout from above.
        display_pretty(f.read(), raw=True)


def run(name):
    get_ipython().system('python {}.py'.format(name))
