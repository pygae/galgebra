"""
Utility Classes
"""

import sys
import collections
import warnings

from io import StringIO  # noqa: F401

# galgebra 0.5.0
warnings.warn(
    "galgebra.utils is deprecated and will be removed. "
    "If you need python 2 compatibility helpers, use a decdicated module like "
    "`six`.", DeprecationWarning, stacklevel=2)

# From https://github.com/benjaminp/six/blob/master/six.py
PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

string_types = str

# https://stackoverflow.com/questions/16176742/python-3-replacement-for-deprecated-compiler-ast-flatten-function


def flatten(x):
    result = []

    for el in x:
        if isinstance(x, collections.Iterable) and not isstr(el):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def isstr(s):
    return isinstance(s, str)
