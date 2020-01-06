"""
Differential operators, for all sympy expressions

For multivector-customized differential operators, see :class:`galgebra.mv.Dop`.
"""
import copy
import numbers
import warnings

from sympy import Symbol, S, Add, simplify, diff, Expr

from . import printer
from . import metric
from . import mv
from .printer import ZERO_STR


def _consolidate_terms(terms):
    """
    Remove zero coefs and consolidate coefs with repeated pdiffs.
    """
    new_coefs = []
    new_pdiffs = []
    for (coef, pd) in terms:
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


################ Scalar Partial Differential Operator Class ############

class Sdop(object):
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

    str_mode = False

    def TSimplify(self):
        return Sdop([
            (metric.Simp.apply(coef), pdiff) for (coef, pdiff) in self.terms
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

    def Sdop_str(self):
        if len(self.terms) == 0:
            return ZERO_STR

        self = self._with_sorted_terms()
        s = ''
        for (coef, pdop) in self.terms:
            coef_str = printer.latex(coef)
            pd_str = printer.latex(pdop)

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

        s = s.replace('+ -','- ')
        s = s[:-3]
        if Sdop.str_mode:
            if len(self.terms) > 1 or isinstance(self.terms[0][0], Add):
                s = '(' + s + ')'
        return s

    def Sdop_latex_str(self):
        if len(self.terms) == 0:
            return ZERO_STR

        self = self._with_sorted_terms()

        s = ''
        for (coef, pdop) in self.terms:
            coef_str = printer.latex(coef)
            pd_str = printer.latex(pdop)
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

        s = s.replace('+ -','- ')
        return s[:-3]

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        return ' ' + latex_str + ' '

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter

        return Printer().doprint(self)

    def __repr__(self):
        return str(self)

    def __init__(self, *args):
        if len(args) == 1 and isinstance(args[0],Symbol):  # Simple Pdop of order 1
            self.terms = ((S(1), Pdop(args[0])),)
        else:
            if len(args) == 2 and isinstance(args[0],list) and isinstance(args[1],list):
                if len(args[0]) != len(args[1]):
                    raise ValueError('In Sdop.__init__ coefficent list and Pdop list must be same length.')
                self.terms = tuple(zip(args[0], args[1]))
            elif len(args) == 1 and isinstance(args[0], (list, tuple)):
                self.terms = tuple(args[0])
            else:
                raise ValueError('In Sdop.__init__ length of args must be 1 or 2 args = '+str(args))

    def __call__(self, arg):
        if isinstance(arg, Sdop):
            terms = []
            for (coef, pdiff) in self.terms:
                new_terms = pdiff(arg.terms)
                new_terms = [(coef * c, p) for c, p in new_terms]
                terms += new_terms
            product = Sdop(terms)
            return Sdop.consolidate_coefs(product)
        else:
            return sum([coef * pdiff(arg) for coef, pdiff in self.terms], S(0))


    def __neg__(self):
        return Sdop([(-coef, pdiff) for coef, pdiff in self.terms])

    @staticmethod
    def Add(sdop1, sdop2):
        if isinstance(sdop1, Sdop) and isinstance(sdop2, Sdop):
            return Sdop(_merge_terms(sdop1.terms, sdop2.terms))
        else:
            # convert values to multiplicative operators
            if isinstance(sdop1, Sdop):
                sdop2 = Sdop([(sdop2, Pdop({}))])
            elif isinstance(sdop2, Sdop):
                sdop1 = Sdop([(sdop1, Pdop({}))])
            else:
                raise TypeError("Neither argument is a Dop instance")
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
        terms = [(sdop * coef, pdiff) for coef, pdiff in self.terms]
        return Sdop(terms)

#################### Partial Derivative Operator Class #################

class Pdop(object):
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

    def __eq__(self,A):
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

    def factor(self):
        """
        If partial derivative operator self.order > 1 factor out first
        order differential operator.  Needed for application of partial
        derivative operator to product of sympy expression and partial
        differential operator.  For example if ``D = Pdop({x:3})`` then::

            (Pdop({x:2}), Pdop({x:1})) = D.factor()
        """
        if self.order == 1:
            return S(0), self
        else:
            new_pdiffs = self.pdiffs.copy()
            x, n = next(iter(new_pdiffs.items()))
            if n == 1:
                del new_pdiffs[x]
            else:
                new_pdiffs[x] -= 1
            return Pdop(new_pdiffs), Pdop(x)

    def __call__(self, arg):
        """
        Calculate nth order partial derivative (order defined by
        self) of :class:`~galgebra.mv.Mv`, :class:`Dop`, :class:`Sdop` or sympy expression
        """
        if self.pdiffs == {}:
            return arg  # result is Pdop identity (1)

        if isinstance(arg, Pdop):  # arg is Pdop
            if arg.pdiffs == {}:  # arg is one
                return self
                #return S(0)  # derivative is zero
            else:  # arg is partial derivative
                pdiffs = copy.copy(arg.pdiffs)
                for key in self.pdiffs:
                    if key in pdiffs:
                        pdiffs[key] += self.pdiffs[key]
                    else:
                        pdiffs[key] = self.pdiffs[key]
            return Pdop(pdiffs)  # result is Pdop

        elif isinstance(arg, mv.Mv):  # arg is multivector
            ga = arg.Ga
            for x in self.pdiffs:
                for i in range(self.pdiffs[x]):
                    arg = ga.pDiff(arg, x)
            return arg  # result is multivector

        elif isinstance(arg, (Expr, Symbol, numbers.Number)):  # arg is sympy expression
            for x in self.pdiffs:
                arg = diff(arg,x,self.pdiffs[x])
            return arg  # derivative is sympy expression

        elif isinstance(arg, (list, tuple)):  # arg is list of tuples (coef, partial derivative)
            terms = list(arg)
            D = self
            while True:
                D, D0 = D.factor()
                for k, term in enumerate(terms):
                    dc = D0(term[0])
                    pd = D0(term[1])
                    #print 'D0, term, dc, pd =', D0, term, dc, pd
                    tmp = []
                    if dc != 0:
                        tmp.append((dc,term[1]))
                    if pd != 0 :
                        tmp.append((term[0],pd))
                    terms[k] = tmp
                terms = [i for o in terms for i in o]  # flatten list one level
                if D == 0:
                    break
            terms = Sdop.consolidate_coefs(terms)
            return terms  # result is list of tuples (coef, partial derivative)
        elif isinstance(arg, Sdop):  # arg is scalar differential operator
            return self(arg.terms)  # result is list of tuples (coef, partial derivative)
        else:
            raise ValueError('In Pdop.__call__ type(arg) = ' + str(type(arg)) + ' not allowed.')

    def __mul__(self, pdop):  # functional product of self and arg (self*arg)
        return self(pdop)

    def __rmul__(self, pdop):  # functional product of arg and self (arg*self)
        if isinstance(pdop, Pdop):
            return pdop(self)
        return Sdop([(pdop, self)])

    def Pdop_str(self):
        if self.order == 0:
            return 'D{}'
        s = 'D'
        for x in self.pdiffs:
            s += '{' + str(x) + '}'
            n = self.pdiffs[x]
            if n > 1:
                s += '^' + str(n)
        return s

    def Pdop_latex_str(self):
        if self.order == 0:
            return ''
        s = r'\frac{\partial'
        if self.order > 1:
            s += '^{' + printer.latex(self.order) + '}'
        s += '}{'
        keys = list(self.pdiffs.keys())
        keys.sort()
        for key in keys:
            i = self.pdiffs[key]
            s += r'\partial ' + printer.latex(key)
            if i > 1:
                s += '^{' + printer.latex(i) + '}'
        s += '}'
        return s

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        return ' ' + latex_str + ' '

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter
        return Printer().doprint(self)

    def __repr__(self):
        return str(self)
