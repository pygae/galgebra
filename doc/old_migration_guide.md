# Migration Guide (Old)

> Note: The APIs have changed since the era of `sympy.galgebra` and `brombo/galgebra`, some properties and methods are deprecated, the supported versions of Python and SymPy have also changed, please check [Changelog](https://galgebra.readthedocs.io/en/latest/changelog.html) and further update your scripts accordingly besides the following. If you encounter any problems, feel free to [open an issue](https://github.com/pygae/galgebra/issues/new)!

# Migrating from [sympy.galgebra](https://docs.sympy.org/0.7.6.1/modules/galgebra/)

GAlgebra is no longer part of SymPy since 1.0.0, if you have an import like this in your source:

```python
from sympy.galgebra.ga import *
```

Simply remove the `sympy.` prefix before `galgebra` then you are good to go:

```python
from galgebra.ga import *
```

# Migrating from [brombo/galgebra](https://github.com/brombo/galgebra)

The `setgapth.py` way to install is now deprecated by `pip install galgebra` and all modules in GAlgebra should be imported from `galgebra`, for example:

```python
from galgebra.printer import Format, Eprint, latex, GaPrinter
from galgebra.ga import Ga
from galgebra.mv import Mv, Nga
```