import re
import string
import sys
import textwrap

import pytest
from IPython.lib.pretty import pretty

from galgebra.ga import lazy_dict, OrderedBiMap
from galgebra._utils import cached_property


class TestOrderedBiMap:

    def test_inverse(self):
        d = OrderedBiMap([('one', 1), ('two', 2)])
        assert d.inverse == OrderedBiMap([(1, 'one'), (2, 'two')])
        assert d.inverse.inverse is d

    def test_annotation(self):
        # check that we can use this like a generic type
        def foo() -> OrderedBiMap[str, int]:
            pass


class TestCachedProperty:

    def test_cache(self):

        class X:
            _val = 0
            @cached_property
            def val(self):
                try:
                    return self._val
                finally:
                    self._val += 1

        # property is computed only when the cache is clear
        x = X()
        assert x.val == 0
        assert x.val == 0
        del x.val
        assert x.val == 1
        assert x.val == 1
        del x.val
        assert x.val == 2
        assert x.val == 2

    @pytest.mark.skipif(
        sys.version_info < (3, 6),
        reason='__set_name__ was not added until Python 3.6'
    )
    def test_annotation(self):
        class X:
            @cached_property
            def val(self) -> 'int':
                pass

        assert X.__annotations__['val'] == 'int'


class TestLazyDict:

    def test_caching(self):
        seen = []
        def f_value(k):
            seen.append(k)
            return len(seen)

        l = lazy_dict({'a': 0}, f_value)
        assert l['a'] == 0
        assert len(seen) == 0

        # should not call the function twice
        assert l['b'] == 1
        assert l['b'] == 1

        assert l['c'] == 2

    def test_repr(self):
        l = lazy_dict({'1': 1}, int)
        assert repr(l) == "lazy_dict({'1': 1}, f_value=<class 'int'>)"

        # needs to be long enough to wrap
        l = lazy_dict({string.ascii_lowercase: string.ascii_lowercase }, int)
        assert pretty(l) == textwrap.dedent("""\
            lazy_dict({'abcdefghijklmnopqrstuvwxyz': 'abcdefghijklmnopqrstuvwxyz'},
                      f_value=<class 'int'>)""")
