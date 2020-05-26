""" Sympy primitives for representing atoms of ga expressions """

from typing import Union

from sympy import Symbol, AtomicExpr, S, Basic, sympify, MatrixExpr
from sympy import Determinant as _Determinant
from sympy.core import numbers
from sympy.core.function import AppliedUndef, UndefinedFunction
from sympy.printing.pretty.stringpict import prettyForm, stringPict
from sympy.printing.pretty.pretty_symbology import U

from . import printer

__all__ = [
    'BasisVectorSymbol',
    'BasisBladeSymbol',
    'BasisBladeNoWedgeSymbol',
    'BasisBaseSymbol',
    'DotProductSymbol',
]


def _all_same(items):
    return all(x == items[0] for x in items)


class BasisVectorSymbol(Symbol):
    """ A symbol representing a basis vector """
    is_commutative = False

    def _sympystr(self, print_obj):
        return printer.Eprint.Base(print_obj._print_Symbol(self))

    def _latex(self, print_obj):
        try:
            return print_obj._print_Symbol(self, style="bold")
        except TypeError:
            # too old a sympy version for `style=`
            return r"\mathbf{{{}}}".format(print_obj._print_Symbol(self))


class _GradedSymbol(AtomicExpr):
    """ Base class for all graded symbols

    Constructing this from a single symbol returns that symbol itself.
    Constructing from no symbols returns the scalar `S.One`.
    This may change in future.
    """
    # the scalar isn't commutative, but __new__ ensures we do not ever create
    # this type of objects for scalars
    is_commutative = False

    def __new__(cls, *args: BasisVectorSymbol) -> Union[
        numbers.One,
        BasisVectorSymbol,
        "_GradedSymbol"
    ]:
        if len(args) == 0:
            return S.One
        elif len(args) == 1:
            return args[0]
        else:
            return super().__new__(cls, *args)


class _JoinedPrinterMixin(Basic):
    """ Helper class to print `Basic.args` joined by symbol.

    Subclasses must populate `_op_sym` and `_op_sym_latex`
    """

    def _sympystr(self, printer):
        return self._op_sympystr.join(
            printer.doprint(v)
            for v in self.args
        )

    def _pretty(self, printer):
        ret = []
        for i, v in enumerate(self.args):
            if i != 0:
                ret.append(self._op_pretty)
            ret.append(printer._print(v))
        return prettyForm(*stringPict.next(*ret))

    def _latex(self, printer):
        return self._op_latex.join(
            printer._print(v)
            for v in self.args
        )


class BasisBaseSymbol(_GradedSymbol, _JoinedPrinterMixin):
    r""" A basis base in a non-orthogonal algebra, such as :math:`e_1 e_2` """
    _op_sympystr = '*'
    _op_pretty = prettyForm('*')
    _op_latex = ''


class BasisBladeSymbol(_GradedSymbol, _JoinedPrinterMixin):
    r""" A basis blade such as :math:`e_1 \wedge e_2` """
    _op_sympystr = '^'
    _op_pretty = prettyForm('^')
    _op_latex = r'\wedge '


class BasisBladeNoWedgeSymbol(BasisBladeSymbol):
    r""" A basis blade with shortened rendering such as :math:`e_{12}` """
    def _split_name(self):
        sub_str = []
        root_str = []
        for basis_vec in self.args:
            split_lst = basis_vec.name.split('_')
            if len(split_lst) != 2:
                raise ValueError('Incompatible basis vector {} for wedgeless printing'.format(basis_vec))
            else:
                sub_str.append(split_lst[1])
                root_str.append(split_lst[0])
        if _all_same(root_str):
            return root_str[0], ''.join(sub_str)
        else:
            raise ValueError('No unique root symbol to use for wedgeless printing')

    def __common_printer(self, printer):
        # print as if we were a basis vector
        root, sub = self._split_name()
        return printer._print(BasisVectorSymbol("{}_{}".format(root, sub)))

    _sympystr = _pretty = _latex = __common_printer


class DotProductSymbol(AtomicExpr):
    """ A symbol used to represent a dot product, like :class:`sympy.DotProduct` """
    is_real = True

    def _sympystr(self, printer):
        a, b = self.args
        return "({}.{})".format(printer._print(a), printer._print(b))

    def _latex(self, printer):
        a, b = self.args
        return r"\left ({}\cdot {}\right ) ".format(printer._print(a), printer._print(b))

    def _pretty(self, printer):
        a, b = self.args
        pform = prettyForm(*stringPict.next(
            printer._print(a),
            printer._print(U('DOT OPERATOR')),
            printer._print(b),
        ))
        return prettyForm(*pform.parens())


class MatrixFunction(UndefinedFunction):
    """ Like a MatrixSymbol, but for functions. """
    def __new__(mcl, name, m, n, **kwargs):
        cls = super().__new__(mcl, name, (AppliedUndef, MatrixExpr), {}, **kwargs)
        cls.shape = sympify(n, strict=True), sympify(n, strict=True)
        return cls


# workaround until sympy/sympy#19354 is merged
if _Determinant.is_commutative is not True:
    class Determinant(_Determinant):
        is_commutative = True
else:
    Determinant = _Determinant
