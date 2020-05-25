class FullSet(set):
    """ A set that contains everything. This is used to trick sympy. """
    def __contains__(self, x):
        return True

    def __iter__(self):
        raise RuntimeError("Set is infinite")

    def __len__(self):
        raise RuntimeError("Set is infinite")

    def __and__(self, other):
        if isinstance(other, set):
            return other
        return NotImplemented
    __rand__ = __and__

    def __or__(self, other):
        if isinstance(other, set):
            return self
        return NotImplemented
    __ror__ = __or__

    def __gt__(self, other):
        if isinstance(other, set):
            return True
        elif isinstance(other, FullSet):
            return False
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, FullSet):
            return False
        return NotImplemented

    def __ge__(self, other):
        return not (self < other)

    def __le__(self, other):
        return not (self > other)

    def __eq__(self, other):
        if isinstance(other, set):
            return False
        elif isinstance(other, FullSet):
            return True
        return NotImplemented

    def __bool__(self):
        return True
