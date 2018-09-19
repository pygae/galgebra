# utils.py

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections

# https://stackoverflow.com/questions/16176742/python-3-replacement-for-deprecated-compiler-ast-flatten-function
def flatten(x):
    result = []

    for el in x:
        if isinstance(x, collections.Iterable) and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result