#!/usr/bin/env python

from distutils.core import setup
from distutils.core import Command
from distutils.command.build_scripts import build_scripts
import sys
import subprocess
import os

modules = [
    'galgebra.ga',
    'galgebra.lt',
    'galgebra.metric',
    'galgebra.mv',
    'galgebra.printer',
    'galgebra.setgapath']

class clean(Command):
    """Cleans *.pyc and debian trashs, so you should get the same copy as
    is in the VCS.
    """

    description = "remove build files"
    user_options = [("all", "a", "the same")]

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        import os
        os.system("find . -name '*.pyc' | xargs rm -f")
        #os.system("rm -f python-build-stamp-2.4")
        #os.system("rm -f MANIFEST")
        os.system("rm -rf build")
        os.system("rm -rf dist")
        os.system('rm -rf doc/GA.log doc/GA.out  doc/GA.synctex.gz doc/GA.toc')

cmdclass = {'clean': clean}

"""
setup(name='galgebra',
      version='0.4.0',
      description='Symbolic Geometric Algebra/Calculus modules for sympy',
      author='Alan Bromborsky',
      author_email='abrombo@verizon.net',
      license='BSD',
      packages=['galgebra','examples','tests','doc'],
      long_description=read('README'),
      classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics'])
"""
