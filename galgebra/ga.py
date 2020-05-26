"""
Geometric Algebra (inherits Metric)
"""
import warnings
import operator
import copy
from itertools import combinations
from functools import reduce
from typing import Tuple, TypeVar, Callable, Dict, Sequence, List, Optional, Union
from ._backports.typing import OrderedDict

from sympy import (
    diff, Rational, Symbol, S, Mul, Add, Expr,
    expand, simplify, eye, trigsimp,
    symbols, sqrt, Matrix,
)

from . import printer
from . import metric
from . import mv
from . import dop
from . import lt
from .atoms import (
    BasisBaseSymbol, BasisBladeSymbol, BasisBladeNoWedgeSymbol,
)
from ._utils import cached_property as _cached_property

half = Rational(1, 2)
one = S(1)
zero = S(0)


# Needed to avoid ambiguity with the methods of the same name, when used in
# type annotations.
_mv = mv
_dop = dop
_lt = lt


# template argument for functions which are Expr -> Expr and Mv -> Mv
_MaybeMv = TypeVar('_MaybeMv', Expr, _mv.Mv)


def is_bases_product(w):
    nc_w = w.args_cnc()
    nc = nc_w[1]
    return len(nc) == 2 or len(nc) == 1 and nc[0].is_Pow and nc[0].exp == 2


_K = TypeVar('_K')
_V = TypeVar('_V')


class lazy_dict(Dict[_K, _V]):
    """
    A dictionary that creates missing entries on the fly.

    When the dictionary is indexed and the key used is not one of the existing
    keys, ``self.f_value(key)`` is called to evaluate the key.  The
    result is then added to the dictionary so that ``self.f_value`` is not
    used to evaluate the same key again.

    Parameters
    ----------
    d :
        Arguments to pass on to the :class:`dict` constructor, typically
        a regular dictionary
    f_value : function
        The function to call to generate a value for a given key
    """
    def __init__(self, d, f_value):
        dict.__init__(self, d)
        self.f_value = f_value

    def __missing__(self, key: _K) -> _V:
        value = self.f_value(key)
        self[key] = value
        return value

    def __repr__(self):
        return '{}({}, f_value={!r})'.format(
            type(self).__qualname__, dict.__repr__(self), self.f_value)

    def _repr_pretty_(self, p, cycle):
        # ipython support
        p_open, p_close = type(self).__qualname__ + '(', ')'
        with p.group(len(p_open), p_open, p_close):
            p.type_pprinters[dict](self, p, cycle)
            p.text(',')
            p.breakable()
            p.text('f_value={}'.format(self.f_value))


def update_and_substitute(expr1, expr2, mul_dict):
    """
    Linear expand expr1 and expr2 to get (summation convention)::

        expr1 = coefs1[i] * bases1[i]
        expr2 = coefs2[j] * bases2[j]

    where ``coefs1`` and ``coefs2`` are lists of are commutative expressions and
    ``bases1`` and ``bases2`` are lists of bases for the geometric algebra.

    Then evaluate::

        expr = coefs1[i] * coefs2[j] * mul_dict[bases1[i], bases2[j]]

    where ``mul_dict[bases1[i], bases2[j]]`` contains the appropriate
    product of ``bases1[i]*bases2[j]`` as a linear combination of scalars and
    bases of the geometric algebra.
    """
    coefs1, bases1 = metric.linear_expand(expr1)
    coefs2, bases2 = metric.linear_expand(expr2)
    expr = S(0)
    for coef1, base1 in zip(coefs1, bases1):
        for coef2, base2 in zip(coefs2, bases2):
            expr += coef1 * coef2 * mul_dict[base1, base2]
    return expr


def nc_subs(expr, base_keys, base_values=None):
    """
    See if expr contains nc (non-commutative) keys in base_keys and substitute corresponding
    value in base_values for nc key.  This was written since standard
    sympy subs was very slow in performing this operation for non-commutative
    keys for long lists of keys.
    """
    if base_values is None:
        [base_keys, base_values] = list(zip(*base_keys))

    if expr.is_commutative:
        return expr
    if isinstance(expr, Add):
        args = expr.args
    else:
        args = [expr]
    s = zero
    for term in args:
        if term.is_commutative:
            s += term
        else:
            c, nc = term.args_cnc(split_1=False)
            key = Mul._from_args(nc)
            coef = Mul._from_args(c)
            if key in base_keys:
                base = base_values[base_keys.index(key)]
                s += coef * base
            else:
                s += term
    return s


_T = TypeVar('_T')
_U = TypeVar('_U')


class GradedTuple(Tuple[Tuple[_T, ...], ...]):
    """ A nested tuple grouped by grade.

    ``self[i]`` refers to a the elements associated with grade ``i``.

    .. attribute:: flat

        :type: Tuple[T]

        The elements flattened out in order of grade.
    """
    def __new__(cls, *args, **kwargs):
        # super does not work here in Python 3.5, as Tuple.__new__ is broken
        self = tuple.__new__(cls, *args, **kwargs)
        self.__dict__['flat'] = tuple(x for single_grade in self for x in single_grade)
        return self

    def __setattr__(self, attr, value):
        raise AttributeError("'GradedTuple' object has no attribute {!r}".format(attr))

    def _map(self, func: Callable[[_T], _U]) -> 'GradedTuple[_U]':
        return GradedTuple(
            tuple(
                func(elem)
                for elem in elems
            )
            for elems in self
        )


class OrderedBiMap(OrderedDict[_K, _V]):
    """ A dict with an ``.inverse`` attribute mapping in the other direction """
    def __init__(self, items):
        # set up the inverse mapping, bypassing our __init__
        self.inverse = OrderedBiMap.__new__(type(self))

        # populate both
        super(OrderedBiMap, self).__init__(items)
        super(OrderedBiMap, self.inverse).__init__([(v, k) for k, v in items])

        # and complete the inverse loop
        self.inverse.inverse = self


class ProductFunction:
    def __init__(self, ga):
        self._ga = ga

    def __call__(self, A: Expr, B: Expr) -> Expr:
        """ Perform the multiplication """
        raise NotImplementedError  # pragma: no cover


class BladeProductFunction(ProductFunction):
    """ Base class for implementations of products between blade representations

    .. automethod:: __call__
    """

    def of_basis_blades(self, blade1: Symbol, blade2: Symbol) -> Expr:
        """ Compute the product of two basis blades """
        raise NotImplementedError  # pragma: no cover

    @_cached_property
    def table_dict(self) -> lazy_dict[Tuple[Symbol, Symbol], Expr]:
        """ A cache of the result of :meth:`of_basis_blades` """
        return lazy_dict({}, f_value=lambda b: self.of_basis_blades(*b))

    def __call__(self, A: Expr, B: Expr) -> Expr:
        return update_and_substitute(A, B, self.table_dict)


class _SingleGradeProductFunction(BladeProductFunction):
    r"""
    Base class for all product functions :math:`\circ` which select a single
    grade from the geometric product,
    :math:`A_r \circ B_s = \left<A_rB_s\right>_{f(r, s)}.
    """
    def _result_grade(self, grade1: int, grade2: int) -> Optional[int]:
        """
        Get the grade to select from the geometric product, for a given
        dot product
        """
        raise NotImplementedError

    def _of_basis_blades_ortho(self, blade1: Symbol, blade2: Symbol):
        # dot (|), left (<), and right (>) products
        # dot product for orthogonal basis
        index1 = self._ga.indexes_to_blades_dict.inverse[blade1]
        index2 = self._ga.indexes_to_blades_dict.inverse[blade2]
        index = list(index1 + index2)

        grade = self._result_grade(len(index1), len(index2))
        if grade is None:
            return zero

        n = len(index)
        sgn = S(1)
        result = S(1)
        ordered = False
        while n > grade:
            ordered = True
            i2 = 1
            while i2 < n:
                i1 = i2 - 1
                index1 = index[i1]
                index2 = index[i2]
                if index1 == index2:
                    n -= 2
                    if n < grade:
                        return zero
                    result *= self._ga.g[index1, index1]
                    index = index[:i1] + index[i2 + 1:]
                elif index1 > index2:
                    ordered = False
                    index[i1] = index2
                    index[i2] = index1
                    sgn = -sgn
                    i2 += 1
                else:
                    i2 += 1
            if ordered:
                break
        if n > grade:
            return zero
        else:
            if index == []:
                return sgn * result
            else:
                return sgn * result * self._ga.indexes_to_blades_dict[tuple(index)]

    def _of_basis_blades_non_ortho(self, blade1: Symbol, blade2: Symbol) -> Expr:
        # dot product of basis blades if basis vectors are non-orthogonal
        # inner (|), left (<), and right (>) products of basis blades
        # grades of input blades
        grade1 = self._ga.blades_to_grades_dict[blade1]
        grade2 = self._ga.blades_to_grades_dict[blade2]

        grade = self._result_grade(grade1, grade2)
        if grade is None:
            return zero

        # Need base rep for blades since that is all we can multiply
        base1 = self._ga.blade_expansion_dict[blade1]
        base2 = self._ga.blade_expansion_dict[blade2]

        # geometric product of basis blades
        base12 = self._ga.mul(base1, base2)
        # blade rep of geometric product
        blade12 = self._ga.base_to_blade_rep(base12)
        # decompose geometric product by grades
        grade_dict = self._ga.grade_decomposition(blade12)

        return grade_dict.get(grade, zero)

    def of_basis_blades(self, blade1: Symbol, blade2: Symbol) -> Expr:
        if self._ga.is_ortho:
            return self._of_basis_blades_ortho(blade1, blade2)
        else:
            return self._of_basis_blades_non_ortho(blade1, blade2)

    def __call__(self, A: Expr, B: Expr) -> Expr:
        return update_and_substitute(A, B, self.table_dict)


class _HestenesDotFunction(_SingleGradeProductFunction):
    def _result_grade(self, grade1: int, grade2: int) -> Optional[int]:
        if grade1 == 0 or grade2 == 0:
            return None
        return abs(grade1 - grade2)


class _ScalarProductFunction(_SingleGradeProductFunction):
    def _result_grade(self, grade1: int, grade2: int) -> Optional[int]:
        return 0


class _LeftContractFunction(_SingleGradeProductFunction):
    def _result_grade(self, grade1: int, grade2: int) -> Optional[int]:
        grade = grade2 - grade1
        if grade < 0:
            return None
        return grade


class _RightContractFunction(_SingleGradeProductFunction):
    def _result_grade(self, grade1: int, grade2: int) -> Optional[int]:
        grade = grade1 - grade2
        if grade < 0:
            return None
        return grade


class _WedgeProductFunction(_SingleGradeProductFunction):
    def _result_grade(self, grade1: int, grade2: int) -> Optional[int]:
        grade = grade1 + grade2
        if grade > self._ga.n:
            return None
        return grade

    # override the base class method with a faster approach
    def of_basis_blades(self, blade1: Symbol, blade2: Symbol) -> Expr:
        # outer (^) product of basis blades
        # this method works for both orthogonal and non-orthogonal basis
        index1 = self._ga.indexes_to_blades_dict.inverse[blade1]
        index2 = self._ga.indexes_to_blades_dict.inverse[blade2]
        index12 = list(index1 + index2)

        grade = self._result_grade(len(index1), len(index2))
        if grade is None:
            return zero

        sgn, wedge12 = Ga.blade_reduce(index12)
        if sgn != 0:
            return sgn * self._ga.indexes_to_blades_dict[tuple(wedge12)]
        else:
            return S(0)


class _GeometricProductFunction(BladeProductFunction):
    def of_basis_blades(self, blade1: Symbol, blade2: Symbol) -> Expr:
        # geometric (*) product for orthogonal basis
        if self._ga.is_ortho:
            index1 = self._ga.indexes_to_blades_dict.inverse[blade1]
            index2 = self._ga.indexes_to_blades_dict.inverse[blade2]
            blade_index = list(index1 + index2)
            repeats = []
            sgn = 1
            for i in range(1, len(blade_index)):
                save = blade_index[i]
                j = i
                while j > 0 and blade_index[j - 1] > save:
                    sgn = -sgn
                    blade_index[j] = blade_index[j - 1]
                    j -= 1
                blade_index[j] = save
                if blade_index[j] == blade_index[j - 1]:
                    repeats.append(save)
            result = S(sgn)
            for i in repeats:
                blade_index.remove(i)
                blade_index.remove(i)
                result *= self._ga.g[i, i]
            if len(blade_index) > 0:
                result *= self._ga.indexes_to_blades_dict[tuple(blade_index)]
            return result
        else:
            base1 = self._ga.blade_to_base_rep(blade1)
            base2 = self._ga.blade_to_base_rep(blade2)
            base12 = self._ga.basic_mul(base1, base2)
            return self._ga.base_to_blade_rep(base12)


class BaseProductFunction(ProductFunction):
    """ Base class for implementations of products between base blade representations """
    def of_basis_bases(self, base1: Symbol, base2: Symbol) -> Expr:
        """ Compute the product of two basis bases """
        raise NotImplementedError


class _BaseGeometricProductFunction(BaseProductFunction):
    def of_basis_bases(self, base1: Symbol, base2: Symbol) -> Expr:
        # geometric product of bases for non-orthogonal basis vectors
        index = self._ga.indexes_to_bases_dict.inverse[base1] + self._ga.indexes_to_bases_dict.inverse[base2]

        coefs, indexes = self._ga.reduce_basis(index)

        return sum((
            coef * self._ga.indexes_to_bases_dict[tuple(index)]
            for coef, index in zip(coefs, indexes)
        ), S(0))

    @_cached_property
    def table_dict(self) -> OrderedDict[Mul, Expr]:
        return OrderedDict(
            (base1 * base2, self.of_basis_bases(base1, base2))
            for base1 in self._ga.bases.flat
            for base2 in self._ga.bases.flat
        )

    def __call__(self, A: Expr, B: Expr) -> Expr:  # geometric product (*) of base representations
        # only multiplicative operation to assume A and B are in base representation
        AxB = expand(A * B)
        AxB = nc_subs(AxB, self.table_dict.items())
        return expand(AxB)


class Ga(metric.Metric):
    r"""
    The vector space (basis, metric, derivatives of basis vectors) is
    defined by the base class :class:`~galgebra.metric.Metric`.

    The instanciating the class :class:`Ga` constructs the geometric algebra of
    the vector space defined by the metric.

    The construction includes the multivector bases, multiplication
    tables or functions for the geometric (``*``), inner (``|``), outer (``^``)
    products, plus the left (``<``) and right (``>``) contractions.  The
    geometric derivative operator and any required connections for the
    derivative are also calculated.

    Except for the geometric product in the case of a non-orthogonal
    set of basis vectors all products and connections (if needed) are
    calculated when needed and place in dictionaries (lists of tuples)
    to be used when needed.  This greatly speeds up evaluations of
    multivector expressions over previous versions of this code since
    the products of multivector bases and connection are not calculated
    unless they are actually needed in the current calculation.

    Only instantiate the :class:`Ga` class via the :class:`~galgebra.mv.Mv` class or any use
    of enhanced printing (text or latex) will cause the bases and multiplication
    table entries to be incorrectly labeled .

    .. rubric:: Inherited from Metric class

    .. autosummary::

        ~galgebra.metric.Metric.g
        ~galgebra.metric.Metric.g_inv
        ~galgebra.metric.Metric.norm
        ~galgebra.metric.Metric.coords
        ~galgebra.metric.Metric.is_ortho
        ~galgebra.metric.Metric.connect_flg
        ~galgebra.metric.Metric.basis
        ~galgebra.metric.Metric.r_symbols
        ~galgebra.metric.Metric.n
        ~galgebra.metric.Metric.n_range
        ~galgebra.metric.Metric.de

    .. rubric:: Basis, basis bases, and basis blades data structures

    .. autosummary::
        ~galgebra.ga.Ga.indexes
        ~galgebra.ga.Ga.bases
        ~galgebra.ga.Ga.blades
        ~galgebra.ga.Ga.mv_blades
        ~galgebra.ga.Ga.coord_vec
        ~galgebra.ga.Ga.indexes_to_blades_dict
        ~galgebra.ga.Ga.indexes_to_bases_dict

    .. rubric:: Multiplication data structures

    The following properties contain implementations of the operators ``*``,
    ``^``, ``|``, ``<``, and ``>``:

    .. autosummary::

        ~galgebra.ga.Ga.mul
        ~galgebra.ga.Ga.wedge
        ~galgebra.ga.Ga.hestenes_dot
        ~galgebra.ga.Ga.left_contract
        ~galgebra.ga.Ga.right_contract

    While behaving like functions, each of these also has a
    :attr:`BladeProductFunction.table_dict` attribute, which contains a lazy lookup table
    of the products of basis blades.

    For non-orthogonal algebras, there is one additional operation, this one mapping
    bases instead of blades. Unlike the others, the ``table_dict`` attribute is
    pre-computed:

    .. autosummary::

        ~galgebra.ga.Ga.basic_mul

    .. rubric:: Reciprocal basis data structures

    .. autosummary::

        ~galgebra.metric.Metric.r_symbols
        ~galgebra.ga.Ga.r_basis
        ~galgebra.ga.Ga.r_basis_dict
        ~galgebra.ga.Ga.r_basis_mv


    .. rubric:: Derivative data structures

    .. attribute:: de

        Derivatives of basis functions.  Two dimensional list. First entry is differentiating coordinate index.
        Second entry is basis vector index.  Quantities are linear combinations of basis vector symbols.

    .. attribute:: grad

        Geometric derivative operator from left. ``grad*F`` returns multivector
        derivative, ``F*grad`` returns differential operator.

    .. attribute:: rgrad

        Geometric derivative operator from right. ``rgrad*F`` returns differential
        operator, ``F*rgrad`` returns multivector derivative.

    .. Sphinx adds all the other members below this docstring

    .. rubric:: Other members

    .. attribute:: dot_mode

        Controls the behavior of :meth:`dot`

        =======  ======================
        value    ``dot`` aliases
        =======  ======================
        ``'|'``  :meth:`hestenes_dot`
        ``'<'``  :meth:`left_contract`
        ``'>'``  :meth:`right_contract`
        =======  ======================
    """

    dual_mode_value = 'I+'
    dual_mode_lst = ['+I', 'I+', '-I', 'I-', '+Iinv', 'Iinv+', '-Iinv', 'Iinv-']

    presets = {'o3d': 'x,y,z:[1,1,1]:[1,1,0]',
               'cyl3d': 'r,theta,z:[1,r**2,1]:[1,1,0]:norm=True',
               'sph3d': 'r,theta,phi:[1,X[0]**2,X[0]**2*cos(X[1])**2]:[1,1,0]:norm=True',
               'para3d': 'u,v,z:[u**2+v**2,u**2+v**2,1]:[1,1,0]:norm=True'}

    @staticmethod
    def dual_mode(mode='I+'):
        """
        Sets mode of multivector dual function for all geometric algebras
        in users program.

        If Ga.dual_mode(mode) not called the default mode is ``'I+'``.

        =====  ============
        mode   return value
        =====  ============
        +I      I*self
        -I     -I*self
        I+      self*I
        I-     -self*I
        +Iinv   Iinv*self
        -Iinv  -Iinv*self
        Iinv+   self*Iinv
        Iinv-  -self*Iinv
        =====  ============
        """
        if mode not in Ga.dual_mode_lst:
            raise ValueError('mode = ' + mode + ' not allowed for Ga.dual_mode.')

        Ga.dual_mode_value = mode

    @staticmethod
    def com(A, B):
        r"""
        Calculate commutator of multivectors :math:`A` and :math:`B`. Returns :math:`(AB-BA)/2`.

        Additionally, commutator and anti-commutator operators are defined by

        .. math::

            \begin{aligned}
                \texttt{A >> B} \equiv & {\displaystyle\frac{AB - BA}{2}} \\
                \texttt{A << B} \equiv & {\displaystyle\frac{AB + BA}{2}}.
            \end{aligned}
        """
        return half * (A * B - B * A)

    @staticmethod
    def build(*args, **kwargs):
        """
        Static method to instantiate geometric algebra and return geometric
        algebra, basis vectors, and grad operator as a tuple.
        """
        GA = Ga(*args, **kwargs)
        basis = list(GA.mv())
        return tuple([GA] + basis)

    @staticmethod
    def preset(setting, root='e', debug=False):

        if setting not in Ga.presets:
            raise ValueError(str(setting) + 'not in Ga.presets.')
        set_lst = Ga.presets[setting].split(':')
        X = symbols(set_lst[0], real=True)
        g = eval(set_lst[1])
        simps = eval(set_lst[2])
        args = [root]
        kwargs = {'g': g, 'coords': X, 'debug': debug, 'I': True, 'gsym': False}

        if len(set_lst) > 3:
            args_lst = set_lst[-1].split(';')
            for arg in args_lst:
                [name, value] = arg.split('=')
                kwargs[name] = eval(value)

        Ga.set_simp(*simps)
        return Ga(*args, **kwargs)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, ga):
        return self.name == ga.name

    def __init__(self, bases, *, wedge=True, **kwargs):
        """
        Parameters
        ----------
        bases :
            Passed as ``basis`` to ``Metric``.
        wedge :
            Use ``^`` symbol to print basis blades
        **kwargs :
            See :class:`galgebra.metric.Metric`.
        """

        self.wedge_print = wedge

        metric.Metric.__init__(self, bases, **kwargs)

        self.par_coords = None

        if self.debug:
            self._print_basis_and_blade_debug()

        self.dot_mode = '|'

        if self.coords is not None:
            self.coords = list(self.coords)

        self.e = mv.Mv(self.blades.flat[-1], ga=self)  # Pseudo-scalar for geometric algebra

        if self.coords is not None:
            self._update_de_from_rbasis()
            self._build_grads()

        if self.connect_flg:
            self._build_connection()

        # Calculate normalized pseudo scalar (I**2 = +/-1)

        self.sing_flg = False

        if self.e_sq.is_number:
            if self.e_sq == S(0):
                self.sing_flg = True
                print('!!!!If I**2 = 0, I cannot be normalized!!!!')
                # raise ValueError('!!!!If I**2 = 0, I cannot be normalized!!!!')
            if self.e_sq > S(0):
                self.i = self.e/sqrt(self.e_sq)
                self.i_inv = self.i
            else:  # I**2 = -1
                self.i = self.e/sqrt(-self.e_sq)
                self.i_inv = -self.i
        else:
            if self.Isq == '+':  # I**2 = 1
                self.i = self.e/sqrt(self.e_sq)
                self.i_inv = self.i
            else:  # I**2 = -1
                self.i = self.e/sqrt(-self.e_sq)
                self.i_inv = -self.i

        if self.debug:
            print('Exit Ga.__init__()')

        self._agrads = {}  # cache of gradient operator with respect to vector a
        self.dslot = -1  # args slot for dervative, -1 for coordinates

        # mystery state used by the Mlt class
        self._mlt_a = []  # List of dummy vectors for Mlt calculations
        self._mlt_acoefs = []  # List of dummy vectors coefficients
        self._mlt_pdiffs = []  # List of lists dummy vector coefficients

        self._XOX = self.mv('XOX', 'vector')  # cached vector for use in is_versor

    @_cached_property
    def coord_vec(self) -> Expr:
        """
        Linear combination of coordinates and basis vectors.  For
        example in orthogonal 3D :math:`x*e_x+y*e_y+z*e_z`.
        """
        if self.coords is None:
            raise ValueError("Ga with no coords has no coord_vec")
        return sum([coord * base for coord, base in zip(self.coords, self.basis)])

    def _reciprocal_of_basis_blade(self, blade: Symbol) -> Expr:
        r"""
        Compute the reciprocal :math:`e^I` of a basis blade :math:`e_I`.

        This is a blade :math:`e^I` such that :math:`\left<\e^Ie_J\right> = \delta_I^J`
        (:cite:`Hestenes`, eq 3.19).
        """
        index = self.indexes_to_blades_dict.inverse[blade]
        r_blade = reduce(self.wedge, [
            self.r_basis[i] for i in index[::-1]
        ], S.One)
        r_blade = r_blade.simplify()

        # normalize at the end
        if not self.is_ortho:
            # r_basis is already normalized if is_ortho is true
            r_blade /= (self.e_sq**len(index))

        return r_blade

    @_cached_property
    def _reciprocal_blade_dict(self) -> lazy_dict:
        """ A dictionary mapping basis blades to their reciprocal blades. """
        return lazy_dict({}, self._reciprocal_of_basis_blade)

    def make_grad(self, a: Union[_mv.Mv, Sequence[Expr]], cmpflg: bool = False) -> mv.Dop:
        r""" Obtain a gradient operator with respect to the multivector a, :math:`\bm{\nabla}_a`."""
        if not isinstance(a, mv.Mv):
            # This might be needed for Mlt, let's leave it till we're sure.
            # Convert to a multivector.
            a = sum((ai * ei for ai, ei in zip(a, self.mv_basis)), self.mv(S.Zero))

        cache_key = (a, cmpflg)

        if cache_key in self._agrads:
            return self._agrads[cache_key]

        # make the grad and cache it
        grad_a = mv.Dop([
            (self.mv(self._reciprocal_blade_dict[base]), dop.Pdop({coef: 1}))
            for coef, base in metric.linear_expand_terms(a.obj)
        ], ga=self, cmpflg=cmpflg)

        self._agrads[cache_key] = grad_a
        return grad_a

    def __str__(self):
        return self.name

    def E(self) -> mv.Mv:  # Unnoromalized pseudo-scalar
        return self.e

    def I(self) -> mv.Mv:  # Noromalized pseudo-scalar
        return self.i

    @property
    def mv_I(self) -> _mv.Mv:
        # This exists for backwards compatibility. Note this is not `I()`!
        # galgebra 0.4.5
        warnings.warn(
            "`ga.mv_I` is deprecated, use `ga.E()` instead, or perhaps `ga.I()`",
            DeprecationWarning, stacklevel=2)
        # default pseudoscalar
        return self.E()

    @property
    def mv_x(self) -> _mv.Mv:
        # This exists for backwards compatibility.
        # galgebra 0.4.5
        warnings.warn(
            "`ga.mv_x` is deprecated, use `ga.mv(your_name, 'vector')` instead",
            DeprecationWarning, stacklevel=2)
        # testing vectors
        return mv.Mv('XxXx', 'vector', ga=self)

    def X(self):
        # galgebra 0.5.0
        warnings.warn(
            "ga.X() is deprecated, use `ga.coord_vec` instead",
            DeprecationWarning, stacklevel=2)
        return self.coord_vec

    @property
    def Pdiffs(self) -> Dict[Symbol, _dop.Pdop]:
        # galgebra 0.4.5
        warnings.warn(
            "ga.Pdiffs[x] is deprecated, use `Pdop(x)` instead",
            DeprecationWarning, stacklevel=2)
        return {x: dop.Pdop(x) for x in self.coords}

    @property
    def sPds(self) -> Dict[Symbol, _dop.Sdop]:
        # galgebra 0.4.5
        warnings.warn(
            "ga.sPds[x] is deprecated, use `Sdop(x)` instead",
            DeprecationWarning, stacklevel=2)
        return {x: dop.Sdop(x) for x in self.coords}

    @property
    def Pdop_identity(self) -> _dop.Pdop:
        # galgebra 0.4.5
        warnings.warn(
            "ga.Pdop_identity is deprecated, use `Pdop({})` instead",
            DeprecationWarning, stacklevel=2)
        return dop.Pdop({})

    @property
    def blades_lst(self) -> List[Symbol]:
        # galgebra 0.5.0
        warnings.warn(
            "ga.blades_lst is deprecated, use `ga.blades.flat[1:]` instead",
            DeprecationWarning, stacklevel=2)
        return list(self.blades.flat[1:])

    @property
    def bases_lst(self) -> List[Symbol]:
        # galgebra 0.5.0
        warnings.warn(
            "ga.bases_lst is deprecated, use `ga.bases.flat[1:]` instead",
            DeprecationWarning, stacklevel=2)
        return list(self.bases.flat[1:])

    @property
    def indexes_lst(self) -> List[Tuple[int, ...]]:
        # galgebra 0.5.0
        warnings.warn(
            "ga.blades_lst is deprecated, use `ga.indexes.flat[1:]` instead",
            DeprecationWarning, stacklevel=2)
        return list(self.indexes.flat[1:])

    def mv(self, root=None, *args, **kwargs) -> Union[_mv.Mv, Tuple[_mv.Mv, ...]]:
        """
        Instanciate and return a multivector for this, 'self',
        geometric algebra.
        """
        if root is None:  # Return ga basis and compute grad and rgrad
            return self.mv_basis

        # ensure that ga is not already in kwargs
        kwargs = dict(ga=self, **kwargs)

        if not isinstance(root, str):
            return mv.Mv(root, *args, **kwargs)

        if ' ' in root and ' ' not in args[0]:
            root_lst = root.split(' ')
            mv_lst = []
            for root in root_lst:
                mv_lst.append(mv.Mv(root, *args, **kwargs))
            return tuple(mv_lst)

        if ' ' in root and ' ' in args[0]:
            root_lst = root.split(' ')
            mvtype_lst = args[0].split(' ')
            if len(root_lst) != len(mvtype_lst):
                raise ValueError('In Ga.mv() for multiple multivectors and ' +
                                 'multivector types incompatible args ' +
                                 str(root_lst) + ' and ' + str(mvtype_lst))

            mv_lst = []
            for root, mv_type in zip(root_lst, mvtype_lst):
                args_list = list(args)
                args_list[0] = mv_type
                args = tuple(args_list)
                mv_lst.append(mv.Mv(root, *args, **kwargs))
            return tuple(mv_lst)

        return mv.Mv(root, *args, **kwargs)

    def mvr(self, norm: bool = True) -> Tuple[_mv.Mv, ...]:
        r"""
        Returns tumple of reciprocal basis vectors.  If norm=True or
        basis vectors are orthogonal the reciprocal basis is normalized
        in the sense that

        .. math:: {i}\cdot e^{j} = \delta_{i}^{j}.

        If the basis is not orthogonal and norm=False then

        .. math:: e_{i}\cdot e^{j} = I^{2}\delta_{i}^{j}.
        """
        if norm and not self.is_ortho:
            return tuple([self.r_basis_mv[i] / self.e_sq for i in self.n_range])
        else:
            return tuple(self.r_basis_mv)

    def bases_dict(self, prefix: str = None) -> Dict[str, Symbol]:
        '''
        returns a dictionary mapping basis element names to their MultiVector
        instances, optionally for specific grades

        if you are lazy,  you might do this to populate your namespace
        with the variables of a given layout.

        >>> locals().update(ga.bases())
        '''
        if prefix is None:
            prefix = 'e'
        bl = self.blades.flat[1:]  # do not include the scalar, which is not named
        var_names = [prefix+''.join([k for k in str(b) if k.isdigit()]) for b in bl]

        return {key: val for key, val in zip(var_names, bl)}

    def _build_grads(self) -> None:
        if not self.is_ortho:
            r_basis = [x / self.e_sq for x in self.r_basis_mv]
        else:
            r_basis = self.r_basis_mv
        if self.norm:
            r_basis = [x / e_norm for x, e_norm in zip(self.r_basis_mv, self.e_norm)]

        pdx = [dop.Pdop(x) for x in self.coords]

        self.grad = mv.Dop(r_basis, pdx, ga=self)
        self.rgrad = mv.Dop(r_basis, pdx, ga=self, cmpflg=True)

    def grads(self) -> Tuple[_mv.Dop, _mv.Dop]:
        if self.coords is None:
            raise ValueError("Ga must have been initialized with coords to compute grads")
        return self.grad, self.rgrad

    def pdop(self, *args, **kwargs) -> _dop.Pdop:
        """ Shorthand to construct a :class:`~galgebra.dop.Pdop` """
        # galgebra 0.4.5
        warnings.warn(
            "`ga.pdop` is deprecated, use `Pdop()` directly.",
            DeprecationWarning, stacklevel=2)
        return dop.Pdop(*args, **kwargs)

    def dop(self, *args, **kwargs) -> _mv.Dop:
        """ Shorthand to construct a :class:`~galgebra.mv.Dop` for this algebra """
        return mv.Dop(*args, ga=self, **kwargs)

    def sdop(self, *args, **kwargs) -> _dop.Sdop:
        """ Shorthand to construct a :class:`~galgebra.dop.Sdop` """
        # galgebra 0.4.5
        warnings.warn(
            "`ga.sdop` is deprecated, use `Sdop()` directly.",
            DeprecationWarning, stacklevel=2)
        return dop.Sdop(*args, **kwargs)

    @property
    def lt_coords(self) -> List[Expr]:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.lt_coords` is deprecated, use the identical `ga.coords`.",
            DeprecationWarning, stacklevel=2)
        return self.coords

    @property
    def lt_x(self) -> Expr:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.lt_x` is deprecated, use the identical `ga.coord_vec`.",
            DeprecationWarning, stacklevel=2)
        return self.coord_vec

    def lt(self, *args, **kwargs):
        """
        Instanciate and return a linear transformation for this, 'self',
        geometric algebra.
        """
        return lt.Lt(*args, ga=self, **kwargs)

    def sm(self, *args, **kwargs) -> 'Sm':
        """
        Instanciate and return a submanifold for this
        geometric algebra.  See :class:`Sm` for instantiation inputs.
        """
        return Sm(*args, ga=self, **kwargs)

    def parametric(self, coords: List[Expr]) -> None:
        if not isinstance(coords, list):
            raise TypeError('In Ga.parametric coords = ' + str(coords) +
                            ' is not a list.')
        if len(coords) != self.n:
            raise ValueError('In Ga.parametric number of parametric functions' +
                             ' not equal to number of coordinates.')

        self.par_coords = {}

        for coord, par_coord in zip(self.coords, coords):
            self.par_coords[coord] = par_coord

    def basis_vectors(self) -> Tuple[Symbol, ...]:
        return tuple(self.basis)

    def _build_basis_base_symbol(self, base_index: Tuple[int, ...]) -> Symbol:
        """ Build a symbol used for the `base_rep` from the given tuple """
        if not base_index:
            return S(1)
        return BasisBaseSymbol(*(self.basis[i] for i in base_index))

    def _build_basis_blade_symbol(self, base_index: Tuple[int, ...]) -> Symbol:
        """ Build a symbol used for the `blade_rep` from the given tuple """
        if self.wedge_print:
            return BasisBladeSymbol(*(self.basis[i] for i in base_index))
        else:
            return BasisBladeNoWedgeSymbol(*(self.basis[i] for i in base_index))

    def _print_basis_and_blade_debug(self) -> None:
        printer.oprint('indexes', self.indexes, 'list(indexes)', self.indexes.flat,
                       'blades', self.blades, 'list(blades)', self.blades.flat,
                       'indexes_to_blades_dict', self.indexes_to_blades_dict,
                       'blades_to_grades_dict', self.blades_to_grades_dict,
                       'blade_super_scripts', self.blade_super_scripts)
        if not self.is_ortho:
            printer.oprint('bases', self.bases, 'list(bases)', self.bases.flat,
                           'indexes_to_bases_dict', self.indexes_to_bases_dict,
                           'bases_to_grades_dict', self.bases_to_grades_dict)

    @_cached_property
    def indexes(self) -> GradedTuple[Tuple[int, ...]]:
        """ Index tuples of basis blades """
        basis_indexes = tuple(self.n_range)
        return GradedTuple(
            tuple(combinations(basis_indexes, i))
            for i in range(len(basis_indexes) + 1)
        )

    @_cached_property
    def blades(self) -> GradedTuple[Symbol]:
        r""" Basis blades symbols by grade.

        The bases for the multivector (geometric) algebra are formed from
        all combinations of the bases of the vector space, including the empty
        combination which is the scalars.

        Each base is represented as a non-commutative symbol of the form

        .. math:: e_{i_{1}}\wedge e_{i_{2}}\wedge ...\wedge e_{i_{r}}.

        where :math:`0 < i_{1} < i_{2} < ... < i_{r}` and :math:`0 < r \le n` the
        dimension of the vector space and :math:`0 < i_{j} \le n`. The total
        number of all symbols of this form is :math:`2^{n}`.

        These are called the blade basis for the geometric algebra and any
        multivector can be represented by a linears combination of these blades.
        The number of basis vectors that are in the symbol for the blade is call
        the grade of the blade.

        Representing the multivector as a linear combination of blades
        gives a blade decomposition of the multivector.

        There is a linear mapping from :attr:`bases` to blades and blades to
        bases so that one can easily convert from one representation to
        another.
        """
        return self.indexes._map(
            lambda index: self._build_basis_blade_symbol(index))

    @_cached_property
    def indexes_to_blades_dict(self) -> OrderedBiMap[Tuple[int, ...], Symbol]:
        """ Bidirectional map from index tuples (:attr:`indices`) to basis blades (:attr:`blades`) """
        return OrderedBiMap(list(zip(self.indexes.flat, self.blades.flat)))

    @_cached_property
    def blades_to_grades_dict(self) -> Dict[Symbol, int]:
        return {
            blade: igrade
            for igrade, grade in enumerate(self.blades)
            for blade in grade
        }

    @_cached_property
    def bases(self) -> GradedTuple[Symbol]:
        r""" Bases (non-commutative sympy symbols) by grade.

        If the basis vectors are not orthogonal a second set of symbols
        is required in addition to the :attr:`blades`, given by:

        .. math:: e_{i_{1}}e_{i_{2}}...e_{i_{r}}

        where :math:`0 < i_{1} < i_{2} < ... < i_{r}` and :math:`0 < r \le n` the
        dimension of the vector space and :math:`0 < i_{j} \le n`. The total
        number of all symbols of this form is :math:`2^{n}`.
        Any multivector can be represented as a linear combination of these bases.

        For the case of an orthogonal set of basis vectors the bases and blades
        are identical, and so this attribute raises :exc:`ValueError`.
        """
        if self.is_ortho:
            raise ValueError("There is no need for bases in orthogonal algebras")
        return self.indexes._map(
            lambda index: self._build_basis_base_symbol(index))

    @_cached_property
    def indexes_to_bases_dict(self) -> OrderedBiMap[Tuple[int, ...], Symbol]:
        """ Bidirectional map from index tuples (:attr:`indices`) to basis bases (:attr:`bases`) """
        return OrderedBiMap(list(zip(self.indexes.flat, self.bases.flat)))

    @_cached_property
    def bases_to_grades_dict(self) -> Dict[Symbol, int]:
        return {
            blade: igrade
            for igrade, grade in enumerate(self.bases)
            for blade in grade
        }

    @_cached_property
    def basis_super_scripts(self) -> List[str]:
        if self.coords is None:
            base0 = str(self.basis[0])
            if '_' in base0:
                sub_index = base0.index('_')
                return [str(base)[sub_index + 1:] for base in self.basis]
            else:
                return [str(i + 1) for i in self.n_range]
        else:
            return [str(coord) for coord in self.coords]

    @_cached_property
    def blade_super_scripts(self) -> GradedTuple[str]:
        return self.indexes._map(lambda base_index: ''.join(
            self.basis_super_scripts[i] for i in base_index
        ))

    @_cached_property
    def mv_blades(self) -> GradedTuple[_mv.Mv]:
        """ :class:`mv.Mv` instances corresponding to :attr:`blades`. """
        return self.blades._map(lambda blade: mv.Mv(blade, ga=self))

    @_cached_property
    def mv_basis(self) -> Tuple[_mv.Mv, ...]:
        """ :class:`mv.Mv` instances corresponding to :attr:`basis`. """
        return tuple(mv.Mv(obj, ga=self) for obj in self.basis)

    @property
    def indexes_to_bases(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.indexes_to_bases` is deprecated, use `ga.indexes_to_bases_dict.items()`",
            DeprecationWarning, stacklevel=2)
        return self.indexes_to_bases_dict.items()

    @property
    def indexes_to_blades(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.indexes_to_blades` is deprecated, use `ga.indexes_to_blades_dict.items()`",
            DeprecationWarning, stacklevel=2)
        return self.indexes_to_blades_dict.items()

    @property
    def bases_to_indexes(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.bases_to_indexes` is deprecated, use `ga.indexes_to_bases_dict.inverse.items()`",
            DeprecationWarning, stacklevel=2)
        return self.indexes_to_bases_dict.inverse.items()

    @property
    def blades_to_indexes(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.blades_to_indexes` is deprecated, use `ga.indexes_to_blades_dict.inverse.items()`",
            DeprecationWarning, stacklevel=2)
        return self.indexes_to_blades_dict.inverse.items()

    @property
    def bases_to_indexes_dict(self) -> OrderedBiMap:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.bases_to_indexes_dict` is deprecated, use `ga.indexes_to_bases_dict.inverse.`",
            DeprecationWarning, stacklevel=2)
        return self.indexes_to_bases_dict.inverse

    @property
    def blades_to_indexes_dict(self) -> OrderedBiMap:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.blades_to_indexes_dict` is deprecated, use `ga.indexes_to_blades_dict.inverse`",
            DeprecationWarning, stacklevel=2)
        return self.indexes_to_blades_dict.inverse

    @property
    def mul_table_dict(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.mul_table_dict` is deprecated, use `ga.mul.table_dict`",
            DeprecationWarning, stacklevel=2)
        return self.mul.table_dict

    @property
    def wedge_table_dict(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.wedge_table_dict` is deprecated, use `ga.wedge.table_dict`",
            DeprecationWarning, stacklevel=2)
        return self.wedge.table_dict

    @property
    def dot_table_dict(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.dot_table_dict` is deprecated, use `ga.hestenes_dot.table_dict`",
            DeprecationWarning, stacklevel=2)
        return self.hestenes_dot.table_dict

    @property
    def left_contract_table_dict(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.left_contract_table_dict` is deprecated, use `ga.left_contract.table_dict`",
            DeprecationWarning, stacklevel=2)
        return self.left_contract.table_dict

    @property
    def right_contract_table_dict(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.right_contract_table_dict` is deprecated, use `ga.right_contract.table_dict`",
            DeprecationWarning, stacklevel=2)
        return self.right_contract.table_dict

    def _build_connection(self):
        # Partial derivatives of multivector bases multiplied (*,^,|,<,>)
        # on left and right (True and False) by reciprocal basis vectors.
        self.connect = {('*', True): [], ('^', True): [], ('|', True): [],
                        ('<', True): [], ('>', True): [], ('*', False): [],
                        ('^', False): [], ('|', False): [], ('<', False): [],
                        ('>', False): []}
        # Partial derivatives of multivector bases
        self._dbases = {}

    ######## Functions for Calculation products of blades/bases ########

    # ******************* Geometric Product (*) ********************** #

    def geometric_product_basis_blades(self, blade12: Tuple[Symbol, Symbol]) -> Expr:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.geometric_product_basis_blades` is deprecated, use `ga.mul.of_basis_blades`",
            DeprecationWarning, stacklevel=2)
        return self.mul.of_basis_blades(*blade12)

    def reduce_basis(self, blst):
        r"""
        Repetitively applies :meth:`reduce_basis_loop` to `blst`
        product representation until normal form is
        realized for non-orthogonal basis

        If the basis vectors are represented by the non-
        commutative symbols :math:`e_1,...,e_n` then a grade :math:`r` base
        is the geometric product :math:`e_{i_1}e_{i_2}\cdots e_{i_r}` where
        :math:`i_1<i_2<\ldots<i_r` (normal form).  Then in galgebra this
        base is represented by a single indexed non-commutative
        symbol with indexes :math:`[i_1,i_2,\ldots,i_r]`.  The total number
        of these bases in an n-dimensional vector space is :math:`2^n`.

        :meth:`reduce_basis` takes the geometric products of basis vectors that
        are not in normal form (out of order) and reduces them to a sum
        of bases that are in normal form (in order).  It does this by
        recursively applying the geometric algebra formula

        .. math::

            e_ie_j = 2(e_i \cdot e_j) - e_je_i

        where the scalar product :math:`e_i \cdot e_j` is obtained from the metric
        tensor of the vector space.  This also allows one to calculate
        the geometric product of any two bases and grade of the
        geometric algebra, and form the multiplication table.
        """
        blst = list(blst)
        if blst == []:  # blst represents scalar
            blst_coef = [1]
            blst_expand = [[]]
            return blst_coef, blst_expand
        blst_expand = [blst]
        blst_coef = [1]
        blst_flg = [False]
        # reduce untill all blst revise flgs are True
        while not reduce(operator.and_, blst_flg):
            for i in range(len(blst_flg)):
                if not blst_flg[i]:  # keep revising if revise flg is False
                    tmp = Ga.reduce_basis_loop(self.g, blst_expand[i])
                    if isinstance(tmp, bool):
                        blst_flg[i] = tmp  # revision of blst_expand[i] complete
                    elif len(tmp) == 3:  # blst_expand[i] contracted
                        blst_coef[i] = tmp[0] * blst_coef[i]
                        blst_expand[i] = tmp[1]
                        blst_flg[i] = tmp[2]
                    else:  # blst_expand[i] revised
                        blst_coef[i] = -blst_coef[i]
                        # if revision force one more pass in case revision
                        # causes repeated index previous to revised pair of
                        # indexes
                        blst_flg[i] = False
                        blst_expand[i] = tmp[3]
                        blst_coef.append(-blst_coef[i] * tmp[0])
                        blst_expand.append(tmp[1])
                        blst_flg.append(tmp[2])
        new_blst_coef = []
        new_blst_expand = []
        for coef, xpand in zip(blst_coef, blst_expand):
            if xpand in new_blst_expand:
                i = new_blst_expand.index(xpand)
                new_blst_coef[i] += coef
            else:
                new_blst_expand.append(xpand)
                new_blst_coef.append(coef)
        return new_blst_coef, new_blst_expand

    @staticmethod
    def reduce_basis_loop(g, blst):
        r"""
        blst is a list of integers :math:`[i_{1},\ldots,i_{r}]` representing the geometric
        product of r basis vectors :math:`a_{{i_1}}\cdots a_{{i_r}}`. :meth:`reduce_basis_loop`
        searches along the list :math:`[i_{1},\ldots,i_{r}]` untill it finds :math:`i_{j} = i_{j+1}`
        and in this case contracts the list, or if :math:`i_{j} > i_{j+1}` it revises
        the list (:math:`\sim i_{j}` means remove :math:`i_{j}` from the list)

        * Case 1: If :math:`i_{j} = i_{j+1}`, return
          :math:`a_{i_{j}}^2` and
          :math:`[i_{1},\ldots,\sim i_{j},\sim i_{j+1},\ldots,i_{r}]`

        * Case 2: If :math:`i_{j} > i_{j+1}`, return
          :math:`a_{i_{j}}.a_{i_{j+1}}`,
          :math:`[i_{1},\ldots,\sim i_{j},\sim i_{j+1},\ldots,i_{r}]`, and
          :math:`[i_{1},\ldots,i_{j+1},i_{j},\ldots,i_{r}]`

        This is an implementation of the formula

        .. math::

            e_i e_j = 2(e_i \cdot e_j) - e_j e_i

        Where :math:`e_i` and :math:`e_j` are basis vectors.
        """
        nblst = len(blst)  # number of basis vectors
        if nblst <= 1:
            return True  # a scalar or vector is already reduced
        for jstep in range(1, nblst):
            istep = jstep - 1
            if blst[istep] == blst[jstep]:  # basis vectorindex is repeated
                i = blst[istep]  # save basis vector index
                if len(blst) > 2:
                    blst = blst[:istep] + blst[jstep + 1:]  # contract blst
                else:
                    blst = []
                if len(blst) <= 1 or jstep == nblst - 1:
                    blst_flg = True  # revision of blst is complete
                else:
                    blst_flg = False  # more revision needed
                return g[i, i], blst, blst_flg
            if blst[istep] > blst[jstep]:  # blst not in normal order
                blst1 = blst[:istep] + blst[jstep + 1:]  # contract blst
                a1 = 2 * g[blst[jstep], blst[istep]]  # coef of contraction
                blst = blst[:istep] + [blst[jstep]] + [blst[istep]] + blst[jstep + 1:]  # revise blst
                if len(blst1) <= 1:
                    blst1_flg = True  # revision of blst is complete
                else:
                    blst1_flg = False  # more revision needed
                return a1, blst1, blst1_flg, blst

        return True  # revision complete, blst in normal order

    # ****************** Outer/wedge (^) product ********************* #

    @staticmethod
    def blade_reduce(lst: List[int]) -> Tuple[int, Optional[List[int]]]:
        """
        Reduce wedge product of basis vectors to normal order.

        `lst` is a list of indicies of basis vectors.  blade_reduce sorts the list
        and determines if the overall number of exchanges in the list is odd or
        even, returning sign changes (``sgn``) and sorted list.  If any two
        indicies in list are equal (wedge product is zero) ``sgn = 0`` and
        ``lst = None`` are returned.
        """
        sgn = S(1)
        for i in range(1, len(lst)):
            save = lst[i]
            j = i
            while j > 0 and lst[j - 1] > save:
                sgn = -sgn
                lst[j] = lst[j - 1]
                j -= 1
            lst[j] = save
            if lst[j] == lst[j - 1]:
                return S(0), None
        return sgn, lst

    def wedge_product_basis_blades(self, blade12: Tuple[Symbol, Symbol]) -> Expr:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.wedge_product_basis_blades` is deprecated, use `ga.wedge.of_basis_blades`",
            DeprecationWarning, stacklevel=2)
        return self.wedge.of_basis_blades(*blade12)

    # ***** Dot (|) product, reft (<) and right (>) contractions ***** #

    def _dot_product_method(self, mode: str) -> _SingleGradeProductFunction:
        if mode == '|':
            return self.hestenes_dot
        elif mode == '<':
            return self.left_contract
        elif mode == '>':
            return self.right_contract
        else:
            raise ValueError('mode={!r} not allowed'.format(mode))

    def dot_product_basis_blades(self, blade12: Tuple[Symbol, Symbol], mode: str) -> Expr:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.dot_product_basis_blades` is deprecated, use `ga.<which-dot>.of_basis_blades` "
            "where `<which-dot>` is one of `hestenes_dot`, `left_contract`, and `right_contract`",
            DeprecationWarning, stacklevel=2)
        return self._dot_product_method(mode)._of_basis_blades_ortho(*blade12)

    def non_orthogonal_dot_product_basis_blades(self, blade12: Tuple[Symbol, Symbol], mode: str) -> Expr:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.non_orthogonal_dot_product_basis_blades` is deprecated, use `ga.<which-dot>.of_basis_blades` "
            "where `<which-dot>` is one of `hestenes_dot`, `left_contract`, and `right_contract`",
            DeprecationWarning, stacklevel=2)
        return self._dot_product_method(mode)._of_basis_blades_non_ortho(*blade12)

    ############# Non-Orthogonal Tables and Dictionaries ###############

    @property
    def basic_mul_table_dict(self) -> OrderedDict[Mul, Expr]:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.basic_mul_table_dict` is deprecated, use `ga.mul.table_dict`",
            DeprecationWarning, stacklevel=2)
        return self.basic_mul.table_dict

    @property
    def basic_mul_table(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.basic_mul_table` is deprecated, use `ga.basic_mul.table_dict.items()`",
            DeprecationWarning, stacklevel=2)
        return list(self.basic_mul.table_dict.items())

    @property
    def basic_mul_keys(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.basic_mul_keys` is deprecated, use `ga.basic_mul.table_dict.keys()`",
            DeprecationWarning, stacklevel=2)
        return list(self.basic_mul.table_dict.keys())

    @property
    def basic_mul_values(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.basic_mul_values` is deprecated, use `ga.basic_mul.table_dict.values()`",
            DeprecationWarning, stacklevel=2)
        return list(self.basic_mul.table_dict.values())

    def non_orthogonal_bases_products(self, base12: Tuple[Symbol, Symbol]) -> Expr:
        # galgebra 0.5.0
        warnings.warn(
            "`ga.non_orthogonal_bases_products` is deprecated, use `ga.basic_mul.of_basis_bases`",
            DeprecationWarning, stacklevel=2)
        return self.basic_mul.of_basis_bases(*base12)

    @_cached_property
    def blade_expansion_dict(self) -> OrderedDict[Symbol, Expr]:
        """ dictionary expanding blade basis in terms of base basis """

        blade_expansion_dict = OrderedDict()

        for blade, index in zip(self.blades.flat, self.indexes.flat):
            grade = len(index)
            if grade <= 1:
                blade_expansion_dict[blade] = blade
            else:
                a = self.indexes_to_blades_dict[index[:1]]
                A = self.indexes_to_blades_dict[index[1:]]
                Aexpand = blade_expansion_dict[A]
                # Formula for outer (^) product of a vector and grade-r multivector
                # a^A_{r} = (a*A + (-1)^{r}*A*a)/2
                # The folowing evaluation takes the most time for setup it is the due to
                # the substitution required for the multiplications
                a_W_A = half * (self.basic_mul(a, Aexpand) - ((-1) ** grade) * self.basic_mul(Aexpand, a))
                blade_expansion_dict[blade] = expand(a_W_A)

        if self.debug:
            print('blade_expansion_dict =', blade_expansion_dict)

        return blade_expansion_dict

    @_cached_property
    def base_expansion_dict(self) -> OrderedDict[Symbol, Expr]:
        """ dictionary expanding base basis in terms of blade basis """
        base_expansion_dict = OrderedDict()

        for base, blade, index in zip(self.bases.flat, self.blades.flat, self.indexes.flat):
            grade = len(index)
            if grade <= 1:
                base_expansion_dict[base] = base
            else:  # back substitution of tridiagonal system
                tmp = self.blade_expansion_dict[blade]
                tmp = tmp.subs(base, -blade)
                tmp = -tmp.subs(base_expansion_dict)
                base_expansion_dict[base] = expand(tmp)

        if self.debug:
            print('base_expansion_dict =', self.base_expansion_dict)

        return base_expansion_dict

    @property
    def base_expansion(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.base_expansion` is deprecated, use `ga.base_expansion_dict.items()`",
            DeprecationWarning, stacklevel=2)
        return list(self.base_expansion_dict.items())

    @property
    def blade_expansion(self):
        # galgebra 0.5.0
        warnings.warn(
            "`ga.blade_expansion` is deprecated, use `ga.blade_expansion_dict.items()`",
            DeprecationWarning, stacklevel=2)
        return list(self.blade_expansion_dict.items())

    def base_to_blade_rep(self, A):

        if self.is_ortho:
            return A
        else:
            # return expand(A).subs(self.base_expansion_dict)
            return nc_subs(expand(A), self.base_expansion_dict.items())

    def blade_to_base_rep(self, A):

        if self.is_ortho:
            return A
        else:
            # return expand(A).subs(self.blade_expansion_dict)
            return nc_subs(expand(A), self.blade_expansion_dict.items())

    ###### Products (*,^,|,<,>) for multivector representations ########

    @_cached_property
    def basic_mul(self) -> BaseProductFunction:
        r""" The geometic product of objects in base form, :math:`A B` """
        if self.is_ortho:
            raise ValueError("No need for this operation for orthogonal algebras")
        return _BaseGeometricProductFunction(self)

    def Mul(self, A: Expr, B: Expr, mode: str = '*') -> Expr:  # Unifies all products into one function
        if mode == '*':
            return self.mul(A, B)
        elif mode == '^':
            return self.wedge(A, B)
        elif mode == '|':
            return self.hestenes_dot(A, B)
        elif mode == '<':
            return self.left_contract(A, B)
        elif mode == '>':
            return self.right_contract(A, B)
        else:
            raise ValueError('Unknown multiplication operator {!r}', mode)

    @_cached_property
    def mul(self) -> BladeProductFunction:
        r""" The geometric product, :math:`A B` """
        return _GeometricProductFunction(self)

    @_cached_property
    def wedge(self) -> BladeProductFunction:
        r""" The wedge product, :math:`A \wedge B` """
        return _WedgeProductFunction(self)

    @_cached_property
    def hestenes_dot(self) -> BladeProductFunction:
        r""" The hestenes dot product, :math:`A \bullet B` """
        return _HestenesDotFunction(self)

    @_cached_property
    def scalar_product(self) -> BladeProductFunction:
        r""" The scalar product, :math:`A * B` """
        return _ScalarProductFunction(self)

    @_cached_property
    def left_contract(self) -> BladeProductFunction:
        r""" The left contraction, :math:`A \rfloor B`  """
        return _LeftContractFunction(self)

    @_cached_property
    def right_contract(self) -> BladeProductFunction:
        r""" The right contraction, :math:`A \lfloor B` """
        return _RightContractFunction(self)

    def dot(self, A: Expr, B: Expr) -> Expr:
        r"""
        Inner product ``|``, ``<``, or ``>``.

        The :attr:`dot_mode` attribute determines which of these is used.
        """
        # forbid something silly like setting dot_mode to the wedge or geometric
        # product
        if self.dot_mode in '^*':
            raise ValueError('"' + str(self.dot_mode) + '" not a legal mode in dot')
        return self.Mul(A, B, mode=self.dot_mode)

    ######################## Helper Functions ##########################

    def grade_decomposition(self, A: _MaybeMv) -> Dict[int, _MaybeMv]:
        """
        Returns dictionary with grades as keys of grades of A.  For example
        if A is a rotor the dictionary keys would be 0 and 2. For a vector
        the single key would be 1.  Note A can be input as a multivector or
        an multivector object (sympy expression).  If A is a multivector the
        dictionary entries are multivectors.  If A is a sympy expression
        (in this case a linear combination of non-commutative symbols) the
        dictionary entries are sympy expressions.
        """
        if isinstance(A, mv.Mv):
            A.blade_rep()
            A.characterise_Mv()
            Aobj = expand(A.obj)
        else:
            Aobj = A
        coefs, blades = metric.linear_expand(Aobj)
        grade_dict = {}
        for coef, blade in zip(coefs, blades):
            if blade == one:
                if 0 in list(grade_dict.keys()):
                    grade_dict[0] += coef
                else:
                    grade_dict[0] = coef
            else:
                grade = self.blades_to_grades_dict[blade]
                if grade in grade_dict:
                    grade_dict[grade] += coef * blade
                else:
                    grade_dict[grade] = coef * blade
        if isinstance(A, mv.Mv):
            for grade in list(grade_dict.keys()):
                grade_dict[grade] = self.mv(grade_dict[grade])
        return grade_dict

    def split_multivector(self, A: _MaybeMv) -> Tuple[Union[Expr, int], Union[Expr, int]]:
        """
        Split multivector :math:`A` into commutative part :math:`a` and
        non-commutative part :math:`A'` so that :math:`A = a+A'`
        """
        if isinstance(A, mv.Mv):
            return self.split_multivector(A.obj)
        else:
            A = expand(A)
            if isinstance(A, Add):
                a = sum([x for x in A.args if x.is_commutative])
                Ap = sum([x for x in A.args if not x.is_commutative])
                return (a, Ap)
            elif isinstance(A, Symbol):
                if A.is_commutative:
                    return (A, 0)
                else:
                    return (0, A)
            else:
                if A.is_commutative:
                    return (A, 0)
                else:
                    return (0, A)

    def remove_scalar_part(self, A: _MaybeMv) -> Union[Expr, int]:
        """
        Return non-commutative part (sympy object) of ``A.obj``.
        """
        if isinstance(A, mv.Mv):
            return self.remove_scalar_part(A.obj)
        else:
            if isinstance(A, Add):
                A = expand(A)
                return sum([x for x in A.args if not x.is_commutative])
            elif isinstance(A, Symbol):
                if A.is_commutative:
                    return 0
                else:
                    return A
            else:
                if A.is_commutative:
                    return 0
                else:
                    return A

    def scalar_part(self, A: _MaybeMv) -> Union[Expr, int]:

        if isinstance(A, mv.Mv):
            return self.scalar_part(A.obj)
        else:
            A = expand(A)
            if isinstance(A, Add):
                return sum([x for x in A.args if x.is_commutative])
            elif isinstance(A, Symbol):
                if A.is_commutative:
                    return A
                else:
                    return 0
            else:
                if A.is_commutative:
                    return A
                else:
                    return 0

    """
        else:
            if A.is_commutative:
                return A
            else:
                return zero
    """

    def grades(self, A: Expr) -> List[int]:  # Return list of grades present in A
        A = self.base_to_blade_rep(A)
        A = expand(A)
        blades = set()
        if isinstance(A, Add):
            args = A.args
        else:
            args = [A]
        for term in args:
            blade = term.args_cnc()[1]
            l_blade = len(blade)
            if l_blade > 0:
                blades.add(blade[0])
            else:
                blades.add(one)
        return sorted({
            self.blades_to_grades_dict[blade]
            for blade in blades
        })

    def reverse(self, A: Expr) -> Expr:  # Calculates reverse of A (see documentation)
        A = expand(A)
        blades = {}
        if isinstance(A, Add):
            args = A.args
        else:
            if A.is_commutative:
                return A
            else:
                args = [A]
        for term in args:
            if term.is_commutative:
                if 0 in blades:
                    blades[0] += term
                else:
                    blades[0] = term
            else:
                _c, nc = term.args_cnc()
                blade = nc[0]
                grade = self.blades_to_grades_dict[blade]
                if grade in blades:
                    blades[grade] += term
                else:
                    blades[grade] = term
        s = zero
        for grade in blades:
            if (grade * (grade - 1)) / 2 % 2 == 0:
                s += blades[grade]
            else:
                s -= blades[grade]
        return s

    def get_grade(self, A: Expr, r: int) -> Expr:  # Return grade r of A, <A>_{r}
        coefs, bases = metric.linear_expand(A)
        return sum((
            coef * base
            for coef, base in zip(coefs, bases)
            if self.blades_to_grades_dict[base] == r
        ), S(0))

    def even_odd(self, A: Expr, even: bool = True) -> Expr:  # Return even or odd part of A
        A = expand(A)
        if A.is_commutative and even:
            return A
        if isinstance(A, Add):
            args = A.args
        else:
            args = [A]
        s = zero
        for term in args:
            if term.is_commutative:
                if even:
                    s += term
            else:
                c, nc = term.args_cnc(split_1=False)
                blade = nc[0]
                grade = self.blades_to_grades_dict[blade]
                if even and grade % 2 == 0:
                    s += Mul._from_args(c) * blade
                elif not even and grade % 2 == 1:
                    s += Mul._from_args(c) * blade
        return s

    @_cached_property
    def e_sq(self) -> Expr:
        r"""
        If ``self.gsym = True`` then :math:`E_{n}^2` is not evaluated, but is represented
        as :math:`E_{n}^2 = (-1)^{n*(n-1)/2}\operatorname{det}(g)` where
        :math:`\operatorname{det}(g)` the determinant
        of the metric tensor can be general scalar function of the coordinates.
        """
        if self.gsym is not None:
            # Define square of pseudo-scalar in terms of metric tensor
            # determinant
            n = self.n
            return (-1) ** (n*(n - 1)//2) * self.detg
        else:
            return simplify(expand((self.e*self.e).scalar()))

    ##################### Multivector derivatives ######################

    @_cached_property
    def r_basis(self) -> List[Expr]:
        r"""
        Reciprocal basis vectors :math:`e^{j}` as linear combination of basis vector symbols.

        These satisfy

        .. math:: e^{j}\cdot e_{k} = \delta_{k}^{j}

        where :math:`\delta_{k}^{j}` is the kronecker delta.  We use the formula
        from Doran and Lasenby 4.94:

        .. math:: e^{j} = (-1)^{j-1}e_{1} \wedge ...e_{j-1} \wedge e_{j+1} \wedge ... \wedge e_{n}*E_{n}^{-1}

        where :math:`E_{n} = e_{1}\wedge ...\wedge e_{n}`.

        For non-orthogonal basis :math:`e^{j}` is not normalized and must be
        divided by :math:`E_{n}^2` (``self.e_sq``) in any relevant calculations.
        """

        if self.debug:
            print('Enter r_basis.\n')

        if self.is_ortho:
            r_basis = [self.basis[i] / self.g[i, i] for i in self.n_range]
        else:

            duals = list(self.blades[self.n - 1])
            # After reverse, the j-th of them is exactly e_{1}^...e_{j-1}^e_{j+1}^...^e_{n}
            duals.reverse()

            sgn = 1
            r_basis = []
            for dual in duals:
                dual_base_rep = self.blade_to_base_rep(dual)
                # {E_n}^{-1} = \frac{E_n}{{E_n}^{2}}
                # r_basis_j = sgn * duals[j] * E_n so it's not normalized, missing a factor of {E_n}^{-2}
                """
                print('blades list =', self.blades.flat)
                print('debug =', expand(self.base_to_blade_rep(self.mul(sgn * dual_base_rep, self.e.obj))))
                print('collect arg =', expand(self.base_to_blade_rep(self.mul(sgn * dual_base_rep, self.e.obj))))
                """
                r_basis_j = metric.collect(expand(self.base_to_blade_rep(self.mul(sgn * dual_base_rep, self.e.obj))), self.blades.flat)
                r_basis.append(r_basis_j)
                # sgn = (-1)**{j-1}
                sgn = -sgn

        if self.debug:
            printer.oprint('E', self.e, 'E**2', self.e_sq, 'unnormalized reciprocal basis =\n', r_basis)
            print('reciprocal basis test =')
            for ei in self.basis:
                for ej in r_basis:
                    ei_dot_ej = self.hestenes_dot(ei, ej)
                    if ei_dot_ej == zero:
                        print('e_{i}|e_{j} = ' + str(ei_dot_ej))
                    else:
                        print('e_{i}|e_{j} = ' + str(expand(ei_dot_ej / self.e_sq)))

        return r_basis

    def _update_de_from_rbasis(self):
        # Replace reciprocal basis vectors with expansion in terms of
        # basis vectors in derivatives of basis vectors.
        de = self.de
        if de is not None:
            for x_i in self.n_range:
                for jb in self.n_range:
                    if not self.is_ortho:
                        de[x_i][jb] = metric.Simp.apply(de[x_i][jb].subs(self.r_basis_dict) / self.e_sq)
                    else:
                        de[x_i][jb] = metric.Simp.apply(de[x_i][jb].subs(self.r_basis_dict))

    @_cached_property
    def g_inv(self) -> Matrix:
        """ inverse of metric tensor, g^{ij} """
        g_inv = eye(self.n)

        for i in self.n_range:
            rx_i = self.r_symbols[i]
            for j in self.n_range:
                rx_j = self.r_symbols[j]
                if j >= i:
                    g_inv[i, j] = self.hestenes_dot(self.r_basis_dict[rx_i], self.r_basis_dict[rx_j])
                    if not self.is_ortho:
                        g_inv[i, j] /= self.e_sq**2
                else:
                    g_inv[i, j] = g_inv[j, i]

        return simplify(g_inv)

    @_cached_property
    def r_basis_dict(self) -> Dict[Symbol, Expr]:
        """ Dictionary to represent reciprocal basis vectors as expansions in terms of basis vectors.

        ``{reciprocal basis symbol: linear combination of basis symbols, ...}``
        """
        return {
            r_symbol: r_base
            for r_symbol, r_base in zip(self.r_symbols, self.r_basis)
        }

    @_cached_property
    def r_basis_mv(self) -> List[_mv.Mv]:
        """ List of reciprocal basis vectors in terms of basis multivectors. """
        return [mv.Mv(r_base, ga=self) for r_base in self.r_basis]

    def er_blade(self, er, blade, mode='*', left=True):
        r"""
        Product (``*``, ``^``, ``|``, ``<``, ``>``) of reciprocal basis vector
        'er' and basis
        blade 'blade' needed for application of derivatives to
        multivectors.  left is 'True' means 'er' is multiplying 'blade'
        on the left, 'False' is for 'er' multiplying 'blade' on the
        right.  Symbolically for left geometric product:

        .. math:: e^{j}*(e_{i_{1}}\wedge ...\wedge e_{i_{r}})
        """
        if mode == '*':
            base = self.blade_to_base_rep(blade)
            if left:
                return self.base_to_blade_rep(self.mul(er, base))
            else:
                return self.base_to_blade_rep(self.mul(base, er))
        elif mode == '^':
            if left:
                return self.wedge(er, blade)
            else:
                return self.wedge(blade, er)
        else:
            if left:
                return self.Mul(er, blade, mode=mode)
            else:
                return self.Mul(blade, er, mode=mode)

    def blade_derivation(self, blade: Symbol, ib: Union[int, Symbol]) -> Expr:
        """
        Calculate derivatives of basis blade 'blade' using derivative of
        basis vectors calculated by metric. 'ib' is the index of the
        coordinate the derivation is with respect to or the coordinate
        symbol.  These are requried for the calculation of the geometric
        derivatives in curvilinear coordinates or for more general
        manifolds.

        'blade_derivation' caches the results in a dictionary, ``self._dbases``,
        so that the derivation for a given blade and coordinate is never
        calculated more that once.

        Note that the return value is not a multivector, but linear combination
        of basis blade symbols.
        """

        if isinstance(ib, int):
            coord = self.coords[ib]
        else:
            coord = ib
            ib = self.coords.index(coord)

        key = (coord, blade)
        if key in self._dbases:
            return self._dbases[key]

        index = self.indexes_to_blades_dict.inverse[blade]
        grade = len(index)

        # differentiate each basis vector separately and sum
        db = S.Zero
        for i in range(grade):
            db += reduce(self.wedge, [
                self.indexes_to_blades_dict[index[:i]],
                self.de[ib][index[i]],
                self.indexes_to_blades_dict[index[i + 1:]]
            ])
        self._dbases[key] = db
        return db

    def pDiff(self, A: _mv.Mv, coord: Union[List, Symbol]) -> _mv.Mv:
        """
        Compute partial derivative of multivector function 'A' with
        respect to coordinate 'coord'.
        """

        if isinstance(coord, list):
            # Perform multiple partial differentiation where coord is
            # a list of differentiation orders for each coordinate and
            # the coordinate is determinded by the list index.  If the
            # element in the list is zero no differentiation is to be
            # performed for that coordinate index.

            dA = copy.copy(A)  # Make copy of A

            for i in self.n_range:
                x = self.coords[i]
                xn = coord[i]
                if xn > 0:  # Differentiate with respect to coordinate x
                    for _j in range(xn):  # xn > 1 multiple differentiation
                        dA = self.pDiff(dA, x)

            return dA

        # Simple partial differentiation, once with respect to a single
        # variable, but including case of non-constant basis vectors

        dA = self.mv(expand(diff(A.obj, coord)))

        if self.connect_flg and self.dslot == -1 and not A.is_scalar():  # Basis blades are function of coordinates
            B = self.remove_scalar_part(A)
            if B != zero:
                if isinstance(B, Add):
                    args = B.args
                else:
                    args = [B]
                for term in args:
                    if not term.is_commutative:
                        c, nc = term.args_cnc(split_1=False)
                        x = self.blade_derivation(nc[0], coord)
                        if x != zero:
                            dA += reduce(operator.mul, c, x)

        return dA

    def grad_sqr(self, A, grad_sqr_mode, mode, left):
        r"""
        Calculate :math:`(grad *_{1} grad) *_{2} A` or :math:`A *_{2} (grad *_{1} grad)`
        where ``grad_sqr_mode`` = :math:`*_{1}` = ``*``, ``^``, or ``|`` and
        ``mode`` = :math:`*_{2}` = ``*``, ``^``, or ``|``.
        """
        Sop, Bop = Ga.DopFop[(grad_sqr_mode, mode)]
        print('(Sop, Bop) =', Sop, Bop)

        print('grad_sqr:A =', A)

        s = zero

        if Sop is False and Bop is False:
            return s

        dA_i = []
        for coord_i in self.coords:
            dA_i.append(self.pDiff(A, coord_i))

        print('dA_i =', dA_i)

        if Sop:
            for i in self.n_range:
                coord_i = self.coords[i]
                if self.connect_flg:
                    s += self.grad_sq_scl_connect[coord_i] * dA_i[i]

                for j in self.n_range:
                    d2A_j = self.pDiff(dA_i[i], self.coords[j])
                    s += self.g_inv[i, j] * d2A_j

        if Bop and self.connect_flg:
            for i in self.n_range:
                coord_i = self.coords[i]
                print('mode =', mode)
                print('dA_i[i] =', dA_i[i])
                if left:
                    if mode == '|':
                        s += self.dot(self.grad_sq_mv_connect[coord_i], dA_i[i])
                    if mode == '^':
                        s += self.wedge(self.grad_sq_mv_connect[coord_i], dA_i[i])
                    if mode == '*':
                        s += self.mul(self.grad_sq_mv_connect[coord_i], dA_i[i])
                else:
                    if mode == '|':
                        s += self.dot(dA_i[i], self.grad_sq_mv_connect[coord_i])
                    if mode == '^':
                        s += self.wedge(dA_i[i], self.grad_sq_mv_connect[coord_i])
                    if mode == '*':
                        s += self.mul(dA_i[i], self.grad_sq_mv_connect[coord_i])
        return s

    def connection(self, rbase, key_base, mode, left):
        """
        Compute required multivector connections of the form
        (Einstein summation convention) :math:`e^{j}*(D_{j}e_{i_{1}...i_{r}})`
        and :math:`(D_{j}e_{i_{1}...i_{r}})*e^{j}` where :math:`*` could be
        ``*``, ``^``, ``|``, ``<``, or ``>`` depending upon the mode, and
        :math:`e^{j}` are reciprocal basis vectors.
        """
        mode_key = (mode, left)
        keys = [i for i, j in self.connect[mode_key]]
        if left:
            key = rbase * key_base
        else:
            key = key_base * rbase
        if key not in keys:
            keys.append(key)
            C = zero
            for ib in self.n_range:
                x = self.blade_derivation(key_base, ib)
                if self.norm:
                    x /= self.e_norm[ib]
                C += self.er_blade(self.r_basis[ib], x, mode, left)
            # Update connection dictionaries
            self.connect[mode_key].append((key, C))
        return C

    def ReciprocalFrame(self, basis: Sequence[_mv.Mv], mode: str = 'norm') -> Tuple[_mv.Mv, ...]:
        r"""
        Compute the reciprocal frame :math:`v^i` of a set of vectors :math:`v_i`.

        Parameters
        ----------
        basis :
            The sequence of vectors :math:`v_i` defining the input frame.
        mode :
            * ``"norm"`` -- indicates that the reciprocal vectors should be
              normalized such that their product with the input vectors is 1,
              :math:`v^i \cdot v_j = \delta_{ij}`.
            * ``"append"`` -- indicates that instead of normalizing, the
              normalization coefficient :math:`E^2` should be appended to the returned tuple.
              One can divide by this coefficient to normalize the vectors.
              The returned vectors are such that
              :math:`v^i \cdot v_j = E^2\delta_{ij}`.

            .. deprecated:: 0.5.0
                Arbitrary strings are interpreted as ``"append"``, but in
                future will be an error
        """
        dim = len(basis)

        indexes = tuple(range(dim))
        index = [()]

        for i in indexes[-2:]:
            index.append(tuple(combinations(indexes, i + 1)))

        MFbasis = []

        for igrade in index[-2:]:
            grade = []
            for iblade in igrade:
                blade = self.mv(S(1), 'scalar')
                for ibasis in iblade:
                    blade ^= basis[ibasis]
                blade = blade.trigsimp()
                grade.append(blade)
            MFbasis.append(grade)
        E = MFbasis[-1][0]
        E_sq = trigsimp((E * E).scalar())

        duals = copy.copy(MFbasis[-2])

        duals.reverse()
        sgn = S(1)
        rbasis = []
        for dual in duals:
            recpv = (sgn * dual * E).trigsimp()
            rbasis.append(recpv)
            sgn = -sgn

        if mode == 'norm':
            for i in range(dim):
                rbasis[i] = rbasis[i] / E_sq
        else:
            if mode != 'append':
                # galgebra 0.5.0
                warnings.warn(
                    "Mode {!r} not understood, falling back to {!r} but this "
                    "is deprecated".format(mode, 'append'),
                    DeprecationWarning, stacklevel=2)
            rbasis.append(E_sq)

        return tuple(rbasis)

    def Mlt(self, *args, **kwargs):
        return lt.Mlt(args[0], self, *args[1:], **kwargs)


class Sm(Ga):
    """
    Submanifold is a geometric algebra defined on a submanifold of a
    base geometric algebra defined on a manifold.  The submanifold is
    defined by a mapping from the coordinates of the base manifold to
    the coordinates of the submanifold. The inputs required to define
    the submanifold are:

    Notes
    -----

    The 'Ga' member function
    'sm' can be used to instantiate the submanifold via (o3d is the base
    manifold)::

        coords = u, v = symbols('u, v', real=True)
        sm_example = o3d.sm([sin(u)*cos(v), sin(u)*sin(v), cos(u)], coords)

        eu, ev = sm_example.mv()
        sm_grad = sm_example.grad
    """
    # __u is to emulate a Python 3.8 positional-only argument, with a clearer
    # spelling than `*args`.
    def __init__(self, __u, __coords, *, ga, norm=False, name=None, root='e', debug=False):
        """
        Parameters
        ----------
        u :
            The coordinate map defining the submanifold
            which is a list of functions of coordinates of the base
            manifold in terms of the coordinates of the submanifold.
            for example if the manifold is a unit sphere then -
            ``u = [sin(u)*cos(v), sin(u)*sin(v), cos(u)]``.

            Alternatively, a parametric vector function
            of the basis vectors of the base manifold.  The
            coefficients of the bases are functions of the coordinates
            (``coords``).  In this case we would call the submanifold
            a "vector" manifold and additional characteristics of the
            manifold can be calculated since we have given an explicit
            embedding of the manifold in the base manifold.
        coords :
            The coordinate list for the submanifold, for
            example ``[u, v]``.
        debug :
            True for debug output
        root : str
            Root symbol for basis vectors
        name : str
            Name of submanifold
        norm : bool
            Normalize basis if True
        ga :
            Base Geometric Algebra
        """

        # print '!!!Enter Sm!!!'

        u = __u
        coords = __coords
        if ga is None:
            raise ValueError('Base geometric algebra must be specified for submanifold.')

        g_base = ga.g_raw
        n_base = ga.n
        n_sub = len(coords)

        # Construct names of basis vectors
        """
        basis_str = ''
        for x in coords:
            basis_str += root + '_' + str(x) + ' '
        basis_str = basis_str[:-1]
        """

        # print 'u =', u

        if isinstance(u, mv.Mv):  # Define vector manifold
            self.ebasis = []
            for coord in coords:
                # Partial derivation of vector function to get basis vectors
                self.ebasis.append(u.diff(coord))

            # print 'sm ebasis =', self.ebasis

            self.g = []
            for b1 in self.ebasis:
                # Metric tensor from dot products of basis vectors
                tmp = []
                for b2 in self.ebasis:
                    tmp.append(b1 | b2)
                self.g.append(tmp)

        else:

            if len(u) != n_base:
                raise ValueError('In submanifold dimension of base manifold' +
                                 ' not equal to dimension of mapping.')
            dxdu = []
            for x_i in u:
                tmp = []
                for u_j in coords:
                    tmp.append(diff(x_i, u_j))
                dxdu.append(tmp)

            # print 'dxdu =', dxdu

            sub_pairs = list(zip(ga.coords, u))

            # Construct metric tensor form coordinate maps
            g = eye(n_sub)  # Zero n_sub x n_sub sympy matrix
            n_range = list(range(n_sub))
            for i in n_range:
                for j in n_range:
                    s = zero
                    for k in ga.n_range:
                        for l in ga.n_range:
                            s += dxdu[k][i] * dxdu[l][j] * g_base[k, l].subs(sub_pairs)
                    g[i, j] = trigsimp(s)

        Ga.__init__(self, root, g=g, coords=coords, norm=norm, debug=debug)

        if isinstance(u, mv.Mv):  # Construct additional functions for vector manifold
            # self.r_basis_mv under construction

            pass

        self.ga = ga
        self.u = u

        if debug:
            print('Exit Sm.__init__()')

    def vpds(self) -> Tuple[_mv.Dop, _mv.Dop]:
        if not self.is_ortho:
            r_basis = [x / self.e_sq for x in self.r_basis_mv]
        else:
            r_basis = self.r_basis_mv
        if self.norm:
            r_basis = [x / e_norm for x, e_norm in zip(self.r_basis_mv, self.e_norm)]

        pdx = [self.Pdiffs[x] for x in self.coords]

        self.vpd = mv.Dop(r_basis, pdx, ga=self)
        self.rvpd = mv.Dop(r_basis, pdx, ga=self, cmpflg=True)
        return self.vpd, self.rvpd
