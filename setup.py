#!/usr/bin/env python
from setuptools import setup, find_packages
from os import path

this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, 'galgebra', '_version.py'), encoding='utf-8') as f:
    exec(f.read())

with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='galgebra',
    version=__version__,  # noqa: F821
    description='Symbolic Geometric Algebra/Calculus package for SymPy.',
    author='Alan Bromborsky',
    maintainer='Hugo Hadfield',
    maintainer_email='hadfield.hugo@gmail.com',
    url='https://github.com/pygae/galgebra',
    license='BSD',
    packages=find_packages(),
    package_dir={'galgebra': 'galgebra'},
    install_requires=['sympy'],
    python_requires='>=3.5.*',
    long_description=long_description,
    long_description_content_type='text/markdown',
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
        'Topic :: Scientific/Engineering :: Physics',
    ],
    project_urls={
        'Documentation': 'http://galgebra.readthedocs.io',
        'Bug Tracker': 'https://github.com/pygae/galgebra/issues',
        'Source Code': 'https://github.com/pygae/galgebra',
    },
)
