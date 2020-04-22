#!/usr/bin/env python
from setuptools import setup, find_packages
import os

version_path = os.path.join('galgebra', '_version.py')
exec(open(version_path).read())

LONG_DESCRIPTION = """
Symbolic Geometric Algebra/Calculus package for SymPy. BSD License.
"""

setup(name='galgebra',
      version=__version__,  # noqa: F821
      description='Symbolic Geometric Algebra/Calculus package for SymPy.',
      author='Alan Bromborsky',
      author_email='hadfield.hugo@gmail.com',
      url='https://github.com/pygae/galgebra',
      license='BSD',
      packages=find_packages(),
      package_dir={'galgebra': 'galgebra'},
      install_requires=['sympy'],
      python_requires='>=3.5.*',
      long_description=LONG_DESCRIPTION,
      classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics'])
