#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.core import Extension
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--noconflict', dest='noconflict', action='store_true', help='Prefix modules with namespace "galgebra"')
parser.set_defaults(noconflict=False)
(args, _) = parser.parse_known_args()
noconflict = vars(args)['noconflict']
# prevent setup() to see the option
try:
    sys.argv.remove('--noconflict')
except Exception as e:
    pass
ga_namespace = 'galgebra' if noconflict else ''

VERSION = '0.4.1.1'
LONG_DESCRIPTION = """
A symbolic geometric algebra module for python. BSD License.
"""

setup(name='galgebra',
      version=VERSION,
      description='Symbolic Geometric Algebra/Calculus modules for sympy',
      author='Alan Bromborsky',
      author_email='hadfield.hugo@gmail.com',
      url='https://github.com/pygae/galgebra',
      license='BSD',
      packages=find_packages(),
      package_dir={'galgebra':'galgebra'},
      install_requires = [
                'sympy',
                'numpy'],
      long_description=LONG_DESCRIPTION,
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
