import sys

# shim for typing.OrderedDict
if sys.version_info >= (3, 7, 2):
    from typing import OrderedDict
else:
    from typing import TypeVar, Mapping
    import collections

    K = TypeVar('K')
    V = TypeVar('V')

    class OrderedDict(collections.OrderedDict, Mapping[K, V]):
        pass

    del K, V, TypeVar, Mapping, collections
