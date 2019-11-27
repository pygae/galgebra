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
    from io import StringIO
else:
    string_types = basestring
    from StringIO import StringIO

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


class _KwargParser:
    """
    Helper function to emulate Python 3 keyword-only arguments.

    Use as::

        def func(x1, **kwargs):
            kw = KwargParser('func', kwargs)
            a = kw.pop('a')
            b = kw.pop('b', 2)
            kw.reject_remaining()
            ...

    To emulate the Python 3 syntax::

        def func(x1, *, a, b=2):
            ...
    """
    def __init__(self, func_name, kwargs):
        self._func_name = func_name
        self._kwargs = kwargs

    def pop(self, arg_name, *default):
        try:
            return self._kwargs.pop(arg_name, *default)
        except KeyError:
            pass
        raise TypeError(
            '{}() missing required keyword-only argument {!r}'
            .format(self._func_name, arg_name)
        )

    def reject_remaining(self):
        if self._kwargs:
            # match the error message to what Python 3 produces
            bad_arg = next(iter(self._kwargs))
            raise TypeError(
                '{}() got an unexpected keyword argument {!r}'
                .format(self._func_name, bad_arg)
            )
