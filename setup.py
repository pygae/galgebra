#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.core import Extension
import os

version_path = os.path.join('galgebra', '_version.py')
exec(open(version_path).read())

LONG_DESCRIPTION = """
Symbolic Geometric Algebra/Calculus package for SymPy. BSD License.
"""

setup(name='galgebra',
      version=__version__,
      description='Symbolic Geometric Algebra/Calculus package for SymPy.',
      author='Alan Bromborsky',
      author_email='hadfield.hugo@gmail.com',
      url='https://github.com/pygae/galgebra',
      license='BSD',
      packages=find_packages(),
      package_dir={'galgebra':'galgebra'},
      install_requires = ['sympy'],
      # 2.7 or >=3.5
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*',
      long_description=LONG_DESCRIPTION,
      classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics'])
