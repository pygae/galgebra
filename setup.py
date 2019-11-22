#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.core import Extension

VERSION = '0.4.1.1'
LONG_DESCRIPTION = """
Symbolic Geometric Algebra/Calculus package for SymPy. BSD License.
"""

setup(name='galgebra',
      version=VERSION,
      description='Symbolic Geometric Algebra/Calculus package for SymPy.',
      author='Alan Bromborsky',
      author_email='hadfield.hugo@gmail.com',
      url='https://github.com/pygae/galgebra',
      license='BSD',
      packages=find_packages(),
      package_dir={'galgebra':'galgebra'},
      install_requires = ['sympy'],
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
