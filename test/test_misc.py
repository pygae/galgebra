from galgebra.ga import OrderedBiMap


class TestOrderedBiMap:

    def test_inverse(self):
        d = OrderedBiMap([('one', 1), ('two', 2)])
        assert d.inverse == OrderedBiMap([(1, 'one'), (2, 'two')])
        assert d.inverse.inverse is d

    def test_annotation(self):
        # check that we can use this like a generic type
        def foo() -> OrderedBiMap[str, int]:
            pass
