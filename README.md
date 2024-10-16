GAlgebra
=========================================

Symbolic Geometric Algebra/Calculus package for SymPy.

[![PyPI](https://img.shields.io/pypi/v/galgebra.svg)](https://pypi.org/project/galgebra/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/galgebra.svg)](https://pypi.org/project/galgebra/)
[![Python CI](https://github.com/pygae/galgebra/actions/workflows/ci.yml/badge.svg)](https://github.com/pygae/galgebra/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/galgebra/badge/?version=latest)](https://galgebra.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/113447311.svg)](https://zenodo.org/badge/latestdoi/113447311)

![](https://raw.githubusercontent.com/pygae/galgebra/master/doc/images/n_vector_positive_spherical.svg?sanitize=true)

Development Status
--------------------

![PyPI - Status](https://img.shields.io/pypi/status/galgebra.svg)
![GitHub contributors](https://img.shields.io/github/contributors/pygae/galgebra.svg)
[![Codecov](https://img.shields.io/codecov/c/github/pygae/galgebra.svg)](https://codecov.io/gh/pygae/galgebra)
[![Maintainability](https://api.codeclimate.com/v1/badges/26d1c1b351d32d2b1097/maintainability)](https://codeclimate.com/github/pygae/galgebra/maintainability)

[brombo/galgebra](https://github.com/brombo/galgebra) was originally written by Alan Bromborsky, but was no longer actively maintained, and as of 2019-11-25 no longer exists.

[pygae/galgebra](https://github.com/pygae/galgebra) is a community fork, maintained by [Pythonic Geometric Algebra Enthusiasts](https://github.com/pygae).

The fork supports Python 3, increases test coverage, sets up CI and linters, maintains releases to [PyPI](https://pypi.org/project/galgebra/#history), improves [docs](http://galgebra.readthedocs.io) and has many bug fixes, see [Changelog](https://galgebra.readthedocs.io/en/latest/changelog.html).

> [!IMPORTANT] 
> Readers of Prof. Alan Macdonald's [Linear and Geometric Algebra](http://www.faculty.luther.edu/~macdonal/laga/index.html) and [Vector and Geometric Calculus](http://www.faculty.luther.edu/~macdonal/vagc/index.html), please check out [**Migrating guide for readers of LAGA&VAGC**](#migrating-guide-for-readers-of-lagavagc) below.
> 
> If you are coming from [sympy.galgebra](https://docs.sympy.org/0.7.6.1/modules/galgebra/) or [brombo/galgebra](https://github.com/brombo/galgebra) (unlikely nowadays), please check out section the old [Migration Guide](https://github.com/pygae/galgebra/blob/master/doc/old_migration_guide.md).

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

```math
\begin{split}\begin{aligned}
  A+B &=  \texttt{A+B} \\
  A-B &=  \texttt{A-B} \\
  AB &=  \texttt{A*B} \\
  A \wedge B &=  \mathtt{A \verb!^! B} \\
  A \cdot B &=  \texttt{A|B} \\
  A \rfloor B &=  \mathtt{A \lt B} \\
  A \lfloor B &=  \mathtt{A \gt B} \\
  A/B &=  \texttt{A/B} \\
\end{aligned}\end{split}
```

### Geometric Calculus

- Geometric Derivative
- Submanifolds
- Linear Transformations
- Differential Operators

The various derivatives of a multivector function is accomplished by multiplying the gradient operator vector with the function:

```math
\begin{aligned}
  \nabla F &=  \texttt{grad*F} \\
  F \bar{\nabla} &=  \texttt{F*rgrad} \\
  \nabla {\wedge}F &=  \mathtt{grad \verb!^! F} \\
  F {\wedge}\bar{\nabla} &=  \mathtt{F \verb!^! rgrad} \\
  \nabla \cdot F &=  \texttt{grad|F} \\
  F \cdot \bar{\nabla} &=  \texttt{F|rgrad} \\
  \nabla \rfloor F &=  \mathtt{grad \lt F} \\
  F \rfloor \bar{\nabla} &=  \mathtt{F \lt rgrad} \\
  \nabla \lfloor F &=  \mathtt{grad \gt F} \\
  F \lfloor \bar{\nabla} &= \mathtt{F \gt rgrad}
\end{aligned}
```

```math
\begin{aligned}
  F \nabla &=  \texttt{F*grad} \\
  \bar{\nabla} F &=  \texttt{rgrad*F} \\
  F {\wedge}\nabla &=  \mathtt{F \verb!^! grad} \\
  \bar{\nabla} {\wedge}F &=  \mathtt{rgrad \verb!^! F} \\
  F \cdot \nabla &=  \texttt{F|grad} \\
  \bar{\nabla}\cdot F &=  \texttt{rgrad|F} \\
  F \rfloor \nabla &=  \mathtt{F \lt grad} \\
  \bar{\nabla} \rfloor F &=  \mathtt{rgrad \lt F} \\
  F \lfloor \nabla &=  \mathtt{F \gt grad} \\
  \bar{\nabla} \lfloor F &= \mathtt{rgrad \gt F}
\end{aligned}
```

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

```math
\begin{aligned}   \langle \mathbf{M} \rangle _3 =& M^{txy}  \boldsymbol{e}_{t}\wedge \boldsymbol{e}_{x}\wedge \boldsymbol{e}_{y} \\  &  + M^{txz}  \boldsymbol{e}_{t}\wedge \boldsymbol{e}_{x}\wedge \boldsymbol{e}_{z} \\  &  + M^{tyz}  \boldsymbol{e}_{t}\wedge \boldsymbol{e}_{y}\wedge \boldsymbol{e}_{z} \\  &  + M^{xyz}  \boldsymbol{e}_{x}\wedge \boldsymbol{e}_{y}\wedge \boldsymbol{e}_{z}  \end{aligned}
```

You may also check out more examples [here](https://github.com/pygae/galgebra/blob/master/examples/).

For detailed documentation, please visit https://galgebra.readthedocs.io/ .

<!-- end: getting-started -->
<!-- begin: installation -->

Installing GAlgebra
---------------------

### Prerequisites

- Works on Linux, Windows, Mac OSX
- [Python](https://www.python.org/) >= 3.8
  - 0.5.0 was the last supported release for Python 3.5-3.7
  - 0.4.x was the last supported release series for Python 2.7
- [SymPy](https://www.sympy.org) >= 1.3 
  - Only SymPy 1.12 is tested via CI, see `.github/workflows/ci.yml` for more details
  - 0.5.0 was the last supported release for SymPy 1.7

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
pytest --nbval examples/ipython/ --nbval examples/primer/ test --nbval-current-env --nbval-sanitize-with test/.nbval_sanitize.cfg
```

This could take more than 10 minutes, please be patient.

<!-- end: installation -->
<!-- begin: migration -->

Migration Guide
----------------

### Migrating guide for readers of LAGA&VAGC

Readers of [Linear and Geometric Algebra](http://www.faculty.luther.edu/~macdonal/laga/index.html) and [Vector and Geometric Calculus](http://www.faculty.luther.edu/~macdonal/vagc/index.html) might be guided by [GAlgebra Primer](http://www.faculty.luther.edu/~macdonal/GAlgebraPrimer.pdf) (version November 29, 2022, accessed May, 2024) to download [GAfiles.zip](http://www.faculty.luther.edu/~macdonal/GAfiles.zip) and copy `gprinter.py`, `lt.py`, `mv.py`, and `GAlgebraInit.py`¸ into where GAlgebra is installed.

These steps are NO LONGER NEEDED since GAlgebra 0.6.0 as they are merge into GAlgebra with tests, copying these files will cause conflicts and regressions of fixed bugs. Instead, you may follow the following steps:

```bash
pip uninstall galgebra
pip install git+https://github.com/pygae/galgebra.git
```

GAlgebra will be installed as `0.6.0-dev` as `0.6.0` has not yet been finalized and published to PyPI.

For minor differences to those files, please check out [the change log for GAlgebra 0.6.0](https://galgebra.readthedocs.io/en/latest/changelog.html#0.6.0). Also please note that:

- `GAlgebraInit.py` is renamed to `primer.py` and can be imported like `from galgebra.primer import *` but it's usage is discouraged, although it saves some boilerplate code, this is not part of GAlgebra's maintained API, GAlgebra might remove it in future.
- Some notebooks from the zip are included in GAlgebra in `examples/primer`.

<!-- end: migration -->

<!-- begin: bundled-resources -->

Bundled Resources
------------------

Note that in the [doc/books](https://github.com/pygae/galgebra/blob/master/doc/books/) directory there are:

- `BookGA.pdf` which is a collection of notes on Geometric Algebra and Calculus based of "Geometric Algebra for Physicists" by Doran and Lasenby and on some papers by Lasenby and Hestenes.
- `galgebra.pdf` which is the original main doc of GAlgebra in PDF format, while the math part is still valid, the part describing the installation and usage of GAlgebra is outdated, please read with caution or visit https://galgebra.readthedocs.io/ instead.
- `Macdonald` which contains bundled supplementary materials for [Linear and Geometric Algebra](http://www.faculty.luther.edu/~macdonal/laga/index.html) and [Vector and Geometric Calculus](http://www.faculty.luther.edu/~macdonal/vagc/index.html) by Alan Macdonald, see [here](https://github.com/pygae/galgebra/blob/master/doc/books/Macdonald/) and [here](https://github.com/pygae/galgebra/blob/master/examples/Macdonald/) for more information.
  - Particularly, `GAlgebraPrimer.pdf` is an archived version of [GAlgebra Primer](http://www.faculty.luther.edu/~macdonal/GAlgebraPrimer.pdf) by Alan Macdonald, last updated on November 29, 2022.

<!-- end: bundled-resources -->

Star History
-------------------

<a href="https://star-history.com/#pygae/galgebra&Date">
 <picture>
   <source media="(prefers-color-scheme: dark)" srcset="https://api.star-history.com/svg?repos=pygae/galgebra&type=Date&theme=dark" />
   <source media="(prefers-color-scheme: light)" srcset="https://api.star-history.com/svg?repos=pygae/galgebra&type=Date" />
   <img alt="Star History Chart" src="https://api.star-history.com/svg?repos=pygae/galgebra&type=Date" />
 </picture>
</a>

Contributors
-------------------

<a href="https://github.com/pygae/galgebra/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=pygae/galgebra" />
</a>

Made with [contrib.rocks](https://contrib.rocks).

Citing This Library
-------------------

For citation information, see [our `CITATION.md` file](https://github.com/pygae/galgebra/blob/master/CITATION.md).
