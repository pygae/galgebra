Galgebra
========

Symbolic Geometric Algebra/Calculus package for SymPy.

[![Build Status](https://travis-ci.com/pygae/galgebra.svg?branch=master)](https://travis-ci.com/pygae/galgebra) [![codecov](https://codecov.io/gh/pygae/galgebra/branch/master/graph/badge.svg)](https://codecov.io/gh/pygae/galgebra)

This package was originally written by Alan Bromborsky.

pygae/galgebra is a community fork, maintained by [Pythonic Geometric Algebra Enthusiasts](https://github.com/pygae).

Note that in the doc directory there is BookGA.pdf which is a collection of notes on 
Geometric Algebra and Calculus based of "Geometric Algebra for Physicists" by Doran and 
Lasenby and on some papers by Lasenby and Hestenes.

Installing Galgebra
---------------------

### Dependencies

- [Python](https://www.python.org/) 2.7 & 3
  - tests pass under both Python 2.7 & 3.6
- [SymPy](https://www.sympy.org)

### Installing Galgebra From PyPI

```bash
pip install galgebra
```

### Installing Galgebra From Source

To install from local source code, run from the repository root:

```bash
pip install .
```

### Running tests to verify the installation

Run from the repository root:

```bash
python setup.py test
```

Or, preferably:

```bash
pip install pytest
pytest test
```

TODO
-----

Get the docs on read the docs

Ensure functions are documented

Set up proper release version system

