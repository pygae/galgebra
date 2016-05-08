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

setup(name='galgebra',
      version='0.4.0',
      description='Symbolic Geometric Algebra/Calculus modules for sympy',
      author='Alan Bromborsky',
      author_email='abrombo@verizon.net',
      license='BSD',
      py_modules=['ga','mv','lt','metric','printer'],
      #packages=find_packages(''),
      package_dir={ga_namespace:'galgebra'},
      install_requires = [
                'sympy',
                'numpy'],
      #long_description=read('README'),
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
