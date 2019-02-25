"""
Utility Classes
"""

import sys
import collections

# From https://github.com/benjaminp/six/blob/master/six.py
PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

if PY3:
    string_types = str
else:
    string_types = basestring

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
    return isinstance(s, string_types)
