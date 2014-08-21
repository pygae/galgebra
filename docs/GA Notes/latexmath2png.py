#!/usr/bin/python
#from __future__ import with_statement # Until Python 2.6
"""
Converts LaTeX math to png images.
Run latexmath2png.py --help for usage instructions.
"""

"""
Author:
    Kamil Kisiel <kamil@kamilkisiel.net>
    URL: http://www.kamilkisiel.net

Revision History:
    2007/04/20 - Initial version

TODO:
    - Make handling of bad input more graceful?
---

Some ideas borrowed from Kjell Fauske's article at http://fauskes.net/nb/htmleqII/

Licensed under the MIT License:

Copyright (c) 2007 Kamil Kisiel <kamil@kamilkisiel.net>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
"""

import os
import sys
import tempfile
import getopt

# Default packages to use when generating output
default_packages = [
        'amsmath',
        'amsthm',
        'amssymb',
        'bm'
        ]

def __build_preamble(packages):
    preamble = '\documentclass{article}\n'
    for p in packages:
        preamble += "\usepackage{%s}\n" % p
    preamble += "\pagestyle{empty}\n\\begin{document}\n"
    return preamble

def __write_output(infile, outdir, workdir = '.', prefix = '', size = 1):
    try:
        # Generate the DVI file
        latexcmd = 'latex -halt-on-error -output-directory %s %s'\
                % (workdir, infile)
        rc = os.system(latexcmd)
        # Something bad happened, abort
        if rc != 0:
            raise Exception('latex error')

        # Convert the DVI file to PNG's
        dvifile = infile.replace('.tex', '.dvi')
        outprefix = os.path.join(outdir, prefix)
        """
        dvicmd = "dvipng -T tight -x %i -z 9 -bg Transparent "\
                "-o %s%%d.png %s" % (size * 1000, outprefix, dvifile)
        """
        dvicmd = 'dvips %s -o %s%%d.ps' % (dvifile,outprefix)
        print dvicmd
        rc = os.system(dvicmd)
        if rc != 0:
            raise Exception('dvipng error')
    finally:
        # Cleanup temporaries
        basefile = infile.replace('.tex', '')
        tempext = [ '.aux', '.dvi', '.log' ]
        for te in tempext:
            tempfile = basefile + te
            if os.path.exists(tempfile):
                os.remove(tempfile)


def math2png(eqs, outdir, packages = default_packages, prefix = '', size = 1):
    """
    Generate png images from $...$ style math environment equations.

    Parameters:
        eqs         - A list of equations
        outdir      - Output directory for PNG images
        packages    - Optional list of packages to include in the LaTeX preamble
        prefix      - Optional prefix for output files
        size        - Scale factor for output
    """
    try:
        # Set the working directory
        workdir = tempfile.gettempdir()

        # Get a temporary file
        fd, texfile = tempfile.mkstemp('.tex', 'eq', workdir, True)

        # Create the TeX document
        tmpstr = ''
        with os.fdopen(fd, 'w+') as f:
            f.write(__build_preamble(packages))
            tmpstr += __build_preamble(packages)
            for eq in eqs:
                f.write("$%s$\n\\newpage\n" % eq)
                tmpstr += "$%s$\n\\newpage\n" % eq
            f.write('\\end{document}')
            tmpstr += '\\end{document}'
            print tmpstr

        __write_output(texfile, outdir, workdir, prefix, size)
    finally:
        if os.path.exists(texfile):
            os.remove(texfile)

def usage():
    print """
Usage: %s [OPTION] ... [FILE] ...
Converts LaTeX math input to PNG.

Options are:
    -h, --help              Display this help information
    --outdir=OUTDIR         PNG file output directory
                            Default: the current working directory
    --packages=PACKAGES     Comma separated list of packages to use
                            Default: amsmath,amsthm,amssymb,bm
    --prefix=PREFIX         Prefix output file names with PREFIX
                            Default: no prefix
    --scale=SCALE           Scale the output by a factor of SCALE.
                            Default: 1 = 100%%

Reads equations from the specified FILEs or standard input if none is given. One
equation is allowed per line of text and each equation is rendered to a separate
PNG image numbered sequentially from 1, with an optional prefix.
    """ % (os.path.split(sys.argv[0])[1])

def main():
    try:
        shortopts = [ 'h', ]
        longopts = [
                'help',
                'outdir=',
                'packages=',
                'prefix=',
                ]
        opts, args = getopt.getopt(sys.argv[1:], shortopts, longopts)
    except getopt.GetoptError, err:
        scriptname = os.path.split(sys.argv[0])[1]
        print "%s: %s" % (scriptname, err)
        print "Try `%s --help` for more information." % scriptname
        sys.exit(2)

    packages = []
    prefix = ''
    outdir = os.getcwd()
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("--packages"):
            packages = a.split(',')
        if o in ("--prefix"):
            prefix = a
        if o in ("--outdir"):
            outdir = os.path.abspath(a)

    input = []
    if args:
        # If filenames were provided on the command line, read their equations
        for a in args:
            fd = os.open(a, os.O_RDONLY)
            with os.fdopen(fd, 'r') as f:
                cur = map(lambda i: i.strip('\n'), f.readlines())
                input.extend(cur)
    else:
        # Otherwise read from stdin
        input = map(lambda i: i.strip('\n'), sys.stdin.readlines())

    # Engage!
    math2png(input, outdir, packages, prefix)

if __name__ == '__main__':
    math2png(['\\Omega = \\beta'],'.',prefix='test',size=5)
    #main()
