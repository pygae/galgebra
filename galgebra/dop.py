"""
Differential operators, for all sympy expressions

For multivector-customized differential operators, see :class:`galgebra.mv.Dop`.
"""
import copy
import numbers
import warnings
from typing import List, Tuple, Any, Iterable

from sympy import Symbol, S, Add, simplify, diff, Expr, Dummy

from . import printer
from . import metric
from .printer import ZERO_STR


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


class Sdop(_BaseDop):
    """
    Scalar differential operator is of the form (Einstein summation)

    .. math:: D = c_{i}*D_{i}

    where the :math:`c_{i}`'s are scalar coefficient (they could be functions)
    and the :math:`D_{i}`'s are partial differential operators (:class:`Pdop`).

    Attributes
    ----------
    terms : tuple of tuple
        the structure :math:`((c_{1},D_{1}),(c_{2},D_{2}), ...)`
    """

    def TSimplify(self):
        return Sdop([
            (metric.Simp.apply(coef), pdiff) for coef, pdiff in self.terms
        ])

    @staticmethod
    def consolidate_coefs(sdop):
        """
        Remove zero coefs and consolidate coefs with repeated pdiffs.
        """
        if isinstance(sdop, Sdop):
            return Sdop(_consolidate_terms(sdop.terms))
        else:
            return _consolidate_terms(sdop)

    def simplify(self, modes=simplify):
        return Sdop([
            (metric.apply_function_list(modes, coef), pdiff)
            for coef, pdiff in self.terms
        ])

    def _with_sorted_terms(self):
        new_terms = sorted(self.terms, key=lambda term: Pdop.sort_key(term[1]))
        return Sdop(new_terms)

    def _sympystr(self, print_obj):
        if len(self.terms) == 0:
            return ZERO_STR

        self = self._with_sorted_terms()
        s = ''
        for coef, pdop in self.terms:
            coef_str = print_obj.doprint(coef)
            pd_str = print_obj.doprint(pdop)

            if coef == S(1):
                s += pd_str
            elif coef == S(-1):
                s += '-' + pd_str
            else:
                if isinstance(coef, Add):
                    s += '(' + coef_str + ')*' + pd_str
                else:
                    s += coef_str + '*' + pd_str
            s += ' + '

        s = s.replace('+ -', '- ')
        s = s[:-3]
        return s

    def _latex(self, print_obj):
        if len(self.terms) == 0:
            return ZERO_STR

        self = self._with_sorted_terms()

        s = ''
        for coef, pdop in self.terms:
            coef_str = print_obj.doprint(coef)
            pd_str = print_obj.doprint(pdop)
            if coef == S(1):
                if pd_str == '':
                    s += '1'
                else:
                    s += pd_str
            elif coef == S(-1):
                if pd_str == '':
                    s += '-1'
                else:
                    s += '-' + pd_str
            else:
                if isinstance(coef, Add):
                    s += r'\left ( ' + coef_str + r'\right ) ' + pd_str
                else:
                    s += coef_str + ' ' + pd_str
            s += ' + '

        s = s.replace('+ -', '- ')
        return s[:-3]

    def __init_from_symbol(self, symbol: Symbol) -> None:
        self.terms = ((S(1), Pdop(symbol)),)

    def __init_from_coef_and_pdiffs(self, coefs: List[Any], pdiffs: List['Pdop']) -> None:
        if not isinstance(coefs, list) or not isinstance(pdiffs, list):
            raise TypeError("coefs and pdiffs must be lists")
        if len(coefs) != len(pdiffs):
            raise ValueError('In Sdop.__init__ coefficent list and Pdop list must be same length.')
        self.terms = tuple(zip(coefs, pdiffs))

    def __init_from_terms(self, terms: Iterable[Tuple[Any, 'Pdop']]) -> None:
        self.terms = tuple(terms)

    def __init__(self, *args):
        if len(args) == 1:
            if isinstance(args[0], Symbol):
                self.__init_from_symbol(*args)
            elif isinstance(args[0], (list, tuple)):
                self.__init_from_terms(*args)
            else:
                raise TypeError(
                    "A symbol or sequence is required (got type {})"
                    .format(type(args[0]).__name__))
        elif len(args) == 2:
            self.__init_from_coef_and_pdiffs(*args)
        else:
            raise TypeError(
                "Sdop() takes from 1 to 2 positional arguments but {} were "
                "given".format(len(args)))

    def __call__(self, arg):
        # Ensure that we return the right type even when there are no terms - we
        # do this by adding `0 * d(arg)/d(nonexistant)`, which must be zero, but
        # will be a zero of the right type.
        dummy_var = Dummy('nonexistant')
        terms = self.terms or ((S(0), Pdop(dummy_var)),)
        return sum([coef * pdiff(arg) for coef, pdiff in terms])

    def __neg__(self):
        return Sdop([(-coef, pdiff) for coef, pdiff in self.terms])

    @staticmethod
    def Add(sdop1, sdop2):
        if isinstance(sdop1, Sdop) and isinstance(sdop2, Sdop):
            return Sdop(_merge_terms(sdop1.terms, sdop2.terms))
        else:
            # convert values to multiplicative operators
            if not isinstance(sdop2, _BaseDop):
                sdop2 = Sdop([(sdop2, Pdop({}))])
            elif not isinstance(sdop1, _BaseDop):
                sdop1 = Sdop([(sdop1, Pdop({}))])
            else:
                return NotImplemented
            return Sdop.Add(sdop1, sdop2)

    def __eq__(self, other):
        if isinstance(other, Sdop):
            diff = self - other
            return len(diff.terms) == 0
        else:
            return NotImplemented

    def __add__(self, sdop):
        return Sdop.Add(self, sdop)

    def __radd__(self, sdop):
        return Sdop.Add(sdop, self)

    def __sub__(self, sdop):
        return Sdop.Add(self, -sdop)

    def __rsub__(self, sdop):
        return Sdop.Add(-self, sdop)

    def __mul__(self, sdopr):
        # alias for applying the operator
        return self.__call__(sdopr)

    def __rmul__(self, sdop):
        return Sdop([(sdop * coef, pdiff) for coef, pdiff in self.terms])

    def _eval_derivative_n_times(self, x, n):
        return Sdop(_eval_derivative_n_times_terms(self.terms, x, n))


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


class Pdop(_BaseDop):
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
            # lower order derivatives first
            self.order,
            # sorted by symbol after that, after expansion
            sorted([
                x.sort_key(order)
                for x, k in self.pdiffs.items()
                for i in range(k)
            ])
        )

    def __eq__(self, A):
        if isinstance(A, Pdop) and self.pdiffs == A.pdiffs:
            return True
        else:
            if len(self.pdiffs) == 0 and A == S(1):
                return True
            return False

    def __init__(self, __arg):
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
            self.pdiffs = __arg
        elif isinstance(__arg, Symbol):  # First order derivative with respect to symbol
            self.pdiffs = {__arg: 1}
        else:
            raise TypeError('A dictionary or symbol is required, got {!r}'.format(__arg))

        self.order = sum(self.pdiffs.values())

    def _eval_derivative_n_times(self, x, n) -> 'Pdop':  # pdiff(self)
        # d is partial derivative
        pdiffs = copy.copy(self.pdiffs)
        if x in pdiffs:
            pdiffs[x] += n
        else:
            pdiffs[x] = n
        return Pdop(pdiffs)

    def __call__(self, arg):
        """
        Calculate nth order partial derivative (order defined by
        self) of expression
        """
        for x, n in self.pdiffs.items():
            arg = _basic_diff(arg, x, n)
        return arg

    def __mul__(self, other):  # functional product of self and arg (self*arg)
        return self(other)

    def __rmul__(self, other):  # functional product of arg and self (arg*self)
        assert not isinstance(other, Pdop)
        return Sdop([(other, self)])

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
