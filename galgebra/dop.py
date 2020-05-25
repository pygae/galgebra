"""
Differential operators, for all sympy expressions

For multivector-customized differential operators, see :class:`galgebra.mv.Dop`.
"""
import abc
import copy
import numbers
import warnings
from typing import List, Tuple, Any, Iterable
import functools
import operator

from sympy import Symbol, S, diff, Expr, Basic
import sympy
from sympy.core.decorators import _sympifyit

from . import printer
from . import _full_set


def _consolidate_terms(terms):
    """
    Remove zero coefs and consolidate coefs with repeated pdiffs.
    """
    new_coefs = []
    new_pdiffs = []
    for coef, pd in terms:
        if coef != S(0):
            if pd in new_pdiffs:
                index = new_pdiffs.index(pd)
                new_coefs[index] += coef
            else:
                new_coefs.append(coef)
                new_pdiffs.append(pd)
    return tuple(zip(new_coefs, new_pdiffs))


def _merge_terms(terms1, terms2):
    """ Concatenate and consolidate two sets of already-consolidated terms """
    pdiffs1 = [pdiff for _, pdiff in terms1]
    pdiffs2 = [pdiff for _, pdiff in terms2]

    pdiffs = pdiffs1 + [x for x in pdiffs2 if x not in pdiffs1]
    coefs = len(pdiffs) * [S(0)]

    for coef, pdiff in terms1:
        index = pdiffs.index(pdiff)
        coefs[index] += coef

    for coef, pdiff in terms2:
        index = pdiffs.index(pdiff)
        coefs[index] += coef

    # remove zeros
    return [(coef, pdiff) for coef, pdiff in zip(coefs, pdiffs) if coef != S(0)]


def _eval_derivative_n_times_terms(terms, x, n):
    for i in range(n):
        new_terms = []
        for k, term in enumerate(terms):
            dc = _basic_diff(term[0], x)
            pd = _basic_diff(term[1], x)
            # print 'D0, term, dc, pd =', D0, term, dc, pd
            if dc != 0:
                new_terms.append((dc, term[1]))
            if pd != 0:
                new_terms.append((term[0], pd))
        terms = new_terms
    return _consolidate_terms(terms)


################ Scalar Partial Differential Operator Class ############

class _BaseDop(printer.GaPrintable):
    """ Base class for differential operators - used to avoid accidental promotion """
    pass


class DiffOpExpr(Expr, _BaseDop):
    free_symbols = _full_set.FullSet()
    is_commutative = False
    _op_priority = 50.0

    def _diff_op_apply(self, x):
        raise NotImplementedError

    @_sympifyit('x', NotImplemented)
    def __add__(self, x):
        return DiffOpAdd(self, x)

    @_sympifyit('x', NotImplemented)
    def __radd__(self, x):
        return DiffOpAdd(x, self)

    @_sympifyit('x', NotImplemented)
    def __mul__(self, x):
        return DiffOpMul(self, x)

    @_sympifyit('x', NotImplemented)
    def __rmul__(self, x):
        return DiffOpMul(x, self)

    def __neg__(self):
        return DiffOpMul(S.NegativeOne, self)

    @_sympifyit('x', NotImplemented)
    def __sub__(self, x):
        return DiffOpAdd(self, -x)

    @_sympifyit('x', NotImplemented)
    def __rsub__(self, x):
        return DiffOpAdd(x, -self)

    def __call__(self, x):
        return self._diff_op_apply(x)


def _DiffOpMul_postprocessor(e: sympy.Mul):
    assert isinstance(e, sympy.Mul)
    return DiffOpMul(*e.args)


def _DiffOpAdd_postprocessor(e: sympy.Add):
    assert isinstance(e, sympy.Add)
    return DiffOpAdd(*e.args)


Basic._constructor_postprocessor_mapping[DiffOpExpr] = {
    "Mul": [_DiffOpMul_postprocessor],
    "Add": [_DiffOpAdd_postprocessor],
}


def _diff_op_ify(x):
    if isinstance(x, DiffOpExpr):
        return x
    elif isinstance(x, sympy.Add):
        return DiffOpAdd(*(a for a in x.args))
    elif isinstance(x, sympy.Mul):
        return DiffOpMul(*(a for a in x.args))
    else:
        return x * DiffOpPartial({})


def _diff_op_apply(d, x):
    if not isinstance(d, DiffOpExpr):
        d = d * DiffOpPartial({})
    return _diff_op_ify(d)._diff_op_apply(x)


class Sdop(abc.ABC):
    @classmethod
    def __subclasshook__(cls, c):
        return issubclass(c, DiffOpExpr)

    @classmethod
    def _from_symbol(cls, symbol: Symbol) -> 'DiffOpExpr':
        return Pdop(symbol)

    @classmethod
    def _from_coef_and_pdiffs(cls, coefs: List[Any], pdiffs: List['Pdop']) -> None:
        if not isinstance(coefs, list) or not isinstance(pdiffs, list):
            raise TypeError("coefs and pdiffs must be lists")
        if len(coefs) != len(pdiffs):
            raise ValueError('In Sdop.__init__ coefficent list and Pdop list must be same length.')
        return cls._from_terms(tuple(zip(coefs, pdiffs)))

    @classmethod
    def _from_terms(cls, terms: Iterable[Tuple[Any, 'Pdop']]) -> None:
        return sum((a * b for a, b in terms), DiffOpAdd.identity)

    def __new__(cls, *args):
        if len(args) == 1:
            if isinstance(args[0], Symbol):
                return cls._from_symbol(*args)
            elif isinstance(args[0], (list, tuple)):
                return cls._from_terms(*args)
            else:
                raise TypeError(
                    "A symbol or sequence is required (got type {})"
                    .format(type(args[0]).__name__))
        elif len(args) == 2:
            return cls._from_coef_and_pdiffs(*args)
        else:
            raise TypeError(
                "Sdop() takes from 1 to 2 positional arguments but {} were "
                "given".format(len(args)))


#################### Partial Derivative Operator Class #################


def _basic_diff(f, x, n=1):
    """ Simple wrapper for `diff` that works for our types too """
    if isinstance(f, (Expr, Symbol, numbers.Number)):  # f is sympy expression
        return diff(f, x, n)
    elif hasattr(f, '_eval_derivative_n_times'):
        # one of our types
        return f._eval_derivative_n_times(x, n)
    else:
        raise ValueError('In_basic_diff type(arg) = ' + str(type(f)) + ' not allowed.')


class DiffOpPartial(DiffOpExpr, sympy.AtomicExpr):
    r"""
    Partial derivative operatorp.

    The partial derivatives are of the form

    .. math::
        \partial_{i_{1}...i_{n}} =
            \frac{\partial^{i_{1}+...+i_{n}}}{\partial{x_{1}^{i_{1}}}...\partial{x_{n}^{i_{n}}}}.

    If :math:`i_{j} = 0` then the partial derivative does not contain the
    :math:`x^{i_{j}}` coordinate.

    Attributes
    ----------
    pdiffs : dict
        A dictionary where coordinates are keys and key value are the number of
        times one differentiates with respect to the key.
    order : int
        Total number of differentiations.
        When this is zero (i.e. when :attr:`pdiffs` is ``{}``) then this object
        is the identity operator, and returns its operand unchanged.
    """

    def sort_key(self, order=None):
        return (
            self.class_key(),
            self.order,
            # sorted by symbol after that, after expansion
            tuple(sorted([
                x.sort_key(order)
                for x, k in self.pdiffs.items()
                for i in range(k)
            ]))
        )

    def __new__(cls, __arg):
        """
        The partial differential operator is a partial derivative with
        respect to a set of real symbols (variables).
        """
        # galgebra 0.4.5
        if __arg is None:
            warnings.warn(
                "`Pdop(None)` is deprecated, use `Pdop({})` instead",
                DeprecationWarning, stacklevel=2)
            __arg = {}

        if isinstance(__arg, dict):  # Pdop defined by dictionary
            pdiffs = __arg
        elif isinstance(__arg, Symbol):  # First order derivative with respect to symbol
            pdiffs = {__arg: 1}
        else:
            raise TypeError('A dictionary or symbol is required, got {!r}'.format(__arg))

        self = super().__new__(cls)
        self.pdiffs = pdiffs
        self.order = sum(self.pdiffs.values())
        return self

    def _eval_derivative(self, x):
        return self._eval_derivative_n_times(x, 1)

    def _eval_derivative_n_times(self, x, n) -> 'Pdop':  # pdiff(self)
        # d is partial derivative
        pdiffs = copy.copy(self.pdiffs)
        if x in pdiffs:
            pdiffs[x] += n
        else:
            pdiffs[x] = n
        return Pdop(pdiffs)

    def _diff_op_apply(self, arg):
        """
        Calculate nth order partial derivative (order defined by
        self) of expression
        """
        for x, n in self.pdiffs.items():
            arg = _basic_diff(arg, x, n)
        return arg

    def _sympystr(self, print_obj):
        if self.order == 0:
            return 'D{}'
        s = 'D'
        for x in self.pdiffs:
            s += '{' + print_obj.doprint(x) + '}'
            n = self.pdiffs[x]
            if n > 1:
                s += '^' + print_obj.doprint(n)
        return s

    def _latex(self, print_obj):
        if self.order == 0:
            return ''
        s = r'\frac{\partial'
        if self.order > 1:
            s += '^{' + print_obj.doprint(self.order) + '}'
        s += '}{'
        keys = list(self.pdiffs.keys())
        keys.sort(key=lambda x: x.sort_key())
        for key in keys:
            i = self.pdiffs[key]
            s += r'\partial ' + print_obj.doprint(key)
            if i > 1:
                s += '^{' + print_obj.doprint(i) + '}'
        s += '}'
        return s

    def __srepr__(self):
        return '{}({})'.format(type(self).__name__, self.pdiffs)

    def _hashable_content(self):
        from sympy.utilities import default_sort_key
        sorted_items = sorted(
            self.pdiffs.items(),
            key=lambda t: (default_sort_key(t[0]), t[1])
        )
        return tuple(sympy.Basic(coeff, sympy.S(n)) for coeff, n in sorted_items)


Pdop = DiffOpPartial


class DiffOpZero(DiffOpExpr):
    def _diff_op_apply(self, x):
        return 0 * x


class DiffOpMul(DiffOpExpr, sympy.Mul):
    identity = DiffOpPartial({})

    def __new__(cls, *args, **kwargs):
        if not args:
            return cls.identity
        pre_coeffs = []
        it_args = iter(args)

        pre_coeffs = []
        diff_ops = []
        diff_operands = []

        # extra pre-multiplied coeffs
        for a in it_args:
            if isinstance(a, DiffOpMul):
                pre_coeffs += a.args[:-1]
                diff_ops.append(a.args[-1])
                break
            if isinstance(a, DiffOpExpr):
                diff_ops.append(a)
                break
            pre_coeffs.append(a)

        # extract differential terms
        for a in it_args:
            if not isinstance(a, DiffOpExpr):
                diff_operands.append(a)
                break
            diff_ops.append(a)

        # must be only one operand
        for a in it_args:
            raise TypeError(
                "Must pass at most one operand after the differential operators"
            )

        # avoid `sympy.Mul` so that this works on multivectors
        if pre_coeffs:
            coeff = functools.reduce(operator.mul, pre_coeffs)
        else:
            coeff = S(1)

        d = cls.identity
        for di in diff_ops:
            d = d._diff_op_apply(di)
        if coeff == S(1):
            self = d
        elif coeff == S(0):
            self = DiffOpZero()
        else:
            self = sympy.Basic.__new__(cls, coeff, d, **kwargs)

        if diff_operands:
            return self._diff_op_apply(diff_operands[0])

        return self

    def _diff_op_apply(self, x):
        coeff, *rest = self.args
        for r in rest[::-1]:
            x = r._diff_op_apply(x)
        return coeff * x

    def diff(self, *args, **kwargs):
        return super().diff(*args, simplify=False, **kwargs)

    def _eval_derivative(self, x):
        coeff, diff = self.args
        return sympy.diff(coeff, x) * diff + coeff * sympy.diff(diff, x)


class DiffOpAdd(DiffOpExpr, sympy.Add):
    identity = DiffOpZero()

    @classmethod
    def _from_args(self, args, is_commutative):
        args = [_diff_op_ify(arg) for arg in args]
        return super()._from_args(args, is_commutative)

    def _diff_op_apply(self, x):
        args = self.args
        assert args
        # avoid `sympy.Add` so that this works on multivectors
        return functools.reduce(operator.add, (_diff_op_apply(a, x) for a in args))

    def _eval_derivative_n_times(self, x, n):
        return DiffOpAdd(*(a.diff(x, n) for a in self.args))

    def _eval_derivative(self, x):
        return DiffOpAdd(*(a.diff(x) for a in self.args))

    def _eval_simplify(self, *args, **kwargs):
        return self

    def diff(self, *args, **kwargs):
        return super().diff(*args, simplify=False, **kwargs)


def _as_terms(d):
    d = d.expand()
    if isinstance(d, DiffOpAdd):
        for a in d.args:
            yield from _as_terms(a)
    elif isinstance(d, DiffOpMul):
        coeff, pdiff = d.args
        for c, p in _as_terms(pdiff):
            yield (coeff*c, p)
    elif isinstance(d, DiffOpPartial):
        yield (S(1), d)
    elif isinstance(d, DiffOpZero):
        pass
    else:
        yield (d, DiffOpMul.identity)
