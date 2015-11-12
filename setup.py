#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.core import Extension


setup(name='galgebra',
      version='0.4.0',
      description='Symbolic Geometric Algebra/Calculus modules for sympy',
      author='Alan Bromborsky',
      author_email='abrombo@verizon.net',
      license='BSD',
      packages=find_packages(),
      package_dir={'galgebra':'galgebra'},
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

