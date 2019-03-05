GAlgebra
=========================================

Symbolic Geometric Algebra/Calculus package for SymPy.

[![PyPI](https://img.shields.io/pypi/v/galgebra.svg)](https://pypi.org/project/galgebra/) [![PyPI - Python Version](https://img.shields.io/pypi/pyversions/galgebra.svg)](https://pypi.org/project/galgebra/) [![Build Status](https://travis-ci.com/pygae/galgebra.svg?branch=master)](https://travis-ci.com/pygae/galgebra) [![Documentation Status](https://readthedocs.org/projects/galgebra/badge/?version=latest)](https://galgebra.readthedocs.io/en/latest/?badge=latest) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/fe7642c639a54d909a36c75db6c2fa49)](https://app.codacy.com/app/utensilcandel/galgebra?utm_source=github.com&utm_medium=referral&utm_content=pygae/galgebra&utm_campaign=Badge_Grade_Settings) [![Codecov](https://img.shields.io/codecov/c/github/pygae/galgebra.svg)](https://codecov.io/gh/pygae/galgebra)

This package was originally written by Alan Bromborsky.

pygae/galgebra is a community fork, maintained by [Pythonic Geometric Algebra Enthusiasts](https://github.com/pygae).

Note that in the doc directory there is BookGA.pdf which is a collection of notes on 
Geometric Algebra and Calculus based of "Geometric Algebra for Physicists" by Doran and 
Lasenby and on some papers by Lasenby and Hestenes.

Installing GAlgebra
---------------------

### Dependencies

- [Python](https://www.python.org/) 2.7 or 3
- [SymPy](https://www.sympy.org)

### Installing GAlgebra From PyPI

```bash
pip install galgebra
```

### Installing GAlgebra From Source

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
