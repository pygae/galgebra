GAlgebra
=========================================

Symbolic Geometric Algebra/Calculus package for SymPy.

[![PyPI](https://img.shields.io/pypi/v/galgebra.svg)](https://pypi.org/project/galgebra/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/galgebra.svg)](https://pypi.org/project/galgebra/)
[![Build Status on CircleCI](https://circleci.com/gh/pygae/galgebra.svg?style=shield)](https://circleci.com/gh/pygae/galgebra)
[![Documentation Status](https://readthedocs.org/projects/galgebra/badge/?version=latest)](https://galgebra.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/113447311.svg)](https://zenodo.org/badge/latestdoi/113447311)

![](https://raw.githubusercontent.com/pygae/galgebra/master/doc/images/n_vector_positive_spherical.svg?sanitize=true)

Development Status
--------------------

![PyPI - Status](https://img.shields.io/pypi/status/galgebra.svg) ![GitHub contributors](https://img.shields.io/github/contributors/pygae/galgebra.svg) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/fe7642c639a54d909a36c75db6c2fa49)](https://app.codacy.com/app/utensilcandel/galgebra?utm_source=github.com&utm_medium=referral&utm_content=pygae/galgebra&utm_campaign=Badge_Grade_Settings) [![Codecov](https://img.shields.io/codecov/c/github/pygae/galgebra.svg)](https://codecov.io/gh/pygae/galgebra)

[brombo/galgebra](https://github.com/brombo/galgebra) was originally written by Alan Bromborsky, but was no longer actively maintained, and as of 2019-11-25 no longer exists.

[pygae/galgebra](https://github.com/pygae/galgebra) is a community fork, maintained by [Pythonic Geometric Algebra Enthusiasts](https://github.com/pygae).

The fork supports Python 3, increases test coverage, set up CI and linters, maintains releases to [PyPI](https://pypi.org/project/galgebra/#history), improves [docs](http://galgebra.readthedocs.io) and has many bug fixes, see [Changelog](https://galgebra.readthedocs.io/en/latest/changelog.html).

Features
--------------------

### Geometric Algebra

- Arbitrary Vector Basis and Metric
- Scalar, Vector, Bivector, Multivector, Pseudoscalar, Spinor, Blade
- Basic Geometic Algebra Operations
  - Sum Difference
  - Geometric Product
  - Outer and Inner Products
  - Left and Right Contractions
  - Reverse, Dual, Exponential
  - Commutator
  - Projection, Reflection, Rotation
  - Reciprocal Frames
- Inspecting Base/Blade Representation
- Symbolic Manipulations
  - `expand`, `factor`, `simplify`, `subs`, `trigsimp` etc.

Overloaded Python operators for basic GA operations:

![](https://raw.githubusercontent.com/pygae/galgebra/master/doc/images/basic_op.svg?sanitize=true)

### Geometric Calculus

- Geometric Derivative
- Submanifolds
- Linear Transformations
- Differential Operators

The various derivatives of a multivector function is accomplished by multiplying the gradient operator vector with the function:

![](https://raw.githubusercontent.com/pygae/galgebra/master/doc/images/grad.svg?sanitize=true) ![](https://raw.githubusercontent.com/pygae/galgebra/master/doc/images/grad_cmp.svg?sanitize=true)

Tip: an example for getting `grad` and `rgrad` of a 3-d Euclidean geometric algebra in rectangular coordinates:

```python
from sympy import symbols
from galgebra.ga import Ga

o3d = Ga('e', g=[1,1,1], coords=symbols('x,y,z',real=True))
(grad,rgrad) = o3d.grads()
```

### Printing

- Enhanced Console Printing
- Latex Printing
  - out-of-the-box support for Jupyter Notebook
  - PDF generation and croping support if you have `pdflatex`/`pdfcrop` installed

<!-- Note: These comments are parsed by our sphinx documentation -->

<!-- begin: getting-started -->

Getting Started
---------------------

After installing GAlgebra (see section [Installing GAlgebra](#installing-galgebra) below), in a Jupyter Notebook:

```python
from sympy import symbols
from galgebra.ga import Ga

from galgebra.printer import Format
Format(Fmode = False, Dmode = True)

st4coords = (t,x,y,z) = symbols('t x y z', real=True)
st4 = Ga('e',
         g=[1,-1,-1,-1],
         coords=st4coords)

M = st4.mv('M','mv',f = True)

M.grade(3).Fmt(3,r'\langle \mathbf{M} \rangle _3')
```

You will see:

![](https://raw.githubusercontent.com/pygae/galgebra/master/doc/images/st4_M3.svg?sanitize=true)

You may also check out more examples [here](https://github.com/pygae/galgebra/blob/master/examples/).

For detailed documentation, please visit https://galgebra.readthedocs.io/ .

**NOTE:** If you are coming from [sympy.galgebra](https://docs.sympy.org/0.7.6.1/modules/galgebra/) or [brombo/galgebra](https://github.com/brombo/galgebra), please check out section [Migration Guide](#migration-guide) below.

<!-- end: getting-started -->
<!-- begin: installation -->

Installing GAlgebra
---------------------

### Prerequisites

- Works on Linux, Windows, Mac OSX
- [Python](https://www.python.org/) >=3.5
- [SymPy](https://www.sympy.org) 1.3, 1.4 or preferably 1.5

### Installing GAlgebra From PyPI (Recommended for users)

```bash
pip install galgebra
```

Then you are all set!

### Installing GAlgebra From Source (Recommended for developers)

To install from the latest source code of GAlgebra:

```bash
git clone https://github.com/pygae/galgebra.git
cd galgebra
pip install -e .
```

Note that the optional `-e` argument is used here for a developer install so modifying the source will take effect immediately without the need of reinstallation.

Now you may run tests to verify the installation, run from the root of the repository:

```bash
pip install pytest
pytest test
```

Further, to run the complete test suite including the ones using [nbval](https://github.com/computationalmodelling/nbval), just run:

```bash
pip install nbval
pytest --nbval examples/ipython/ test --current-env --sanitize-with test/.nbval_sanitize.cfg
```

This could take more than 10 minutes, please be patient.

<!-- end: installation -->
<!-- begin: migration -->

Migration Guide
----------------

### Migrating from [sympy.galgebra](https://docs.sympy.org/0.7.6.1/modules/galgebra/)

GAlgebra is no longer part of SymPy since 1.0.0, if you have an import like this in your source:

```python
from sympy.galgebra.ga import *
```

Simply remove the `sympy.` prefix before `galgebra` then you are good to go:

```python
from galgebra.ga import *
```

### Migrating from [brombo/galgebra](https://github.com/pygae/galgebra)

The `setgapth.py` way to install is now deprecated by `pip install galgebra` and all modules in GAlgebra should be imported from `galgebra`, for example:

```python
from galgebra.printer import Format, Eprint, Get_Program, latex, GaPrinter
from galgebra.ga import Ga, one, zero
from galgebra.mv import Mv, Nga
```

<!-- end: migration -->
<!-- begin: bundled-resources -->

Bundled Resources
------------------

Note that in the [doc/books](https://github.com/pygae/galgebra/blob/master/doc/books/) directory there are:

- `BookGA.pdf` which is a collection of notes on Geometric Algebra and Calculus based of "Geometric Algebra for Physicists" by Doran and Lasenby and on some papers by Lasenby and Hestenes.
- `galgebra.pdf` which is the original main doc of GAlgebra in PDF format, while the math part is still valid, the part describing the installation and usage of GAlgebra is outdated, please read with cautious or visit https://galgebra.readthedocs.io/ instead.
- `Macdonald` which constains bundled supplementary materials for [Linear and Geometric Algebra](http://www.faculty.luther.edu/~macdonal/laga/index.html) and [Vector and Geometric Calculus](http://www.faculty.luther.edu/~macdonal/vagc/index.html) by Alan Macdonald, see [here](https://github.com/pygae/galgebra/blob/master/doc/books/Macdonald/) and [here](https://github.com/pygae/galgebra/blob/master/examples/Macdonald/) for more information.

<!-- end: bundled-resources -->
