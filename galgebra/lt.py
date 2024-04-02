"""
Multivector Linear Transformation
"""

import inspect
import types
import itertools
import warnings
from copy import copy
from functools import reduce
from typing import Mapping

from sympy import (
    expand, symbols, Matrix, Transpose, zeros, Symbol, Function, S, Add, Expr, simplify
)
from sympy.printing.latex import LatexPrinter as _LatexPrinter
from sympy.printing.str import StrPrinter as _StrPrinter

from ._utils import cached_property as _cached_property

from . import printer
from . import metric
from . import mv


# Add custom settings to the builtin latex printer
_LatexPrinter._default_settings.update({
    'galgebra_mlt_lcnt': 1
})
_StrPrinter._default_settings.update({
    'galgebra_mlt_lcnt': 1
})


# ## GSG code starts ###
def Symbolic_Matrix(kernel, coords=None, f=False, mode='g'):
    """
    Returns a square real matrix the entries of which are symbolic
    constants or symbolic functions of the coordinates.
    - `kernel` is a one-letter string.  It specifies the kernel letter of
      indexed symbols or functions used to specify the matrix's entries
    - `coords` is a list or tuple.  Its entries are used to label the
      components of a vector.
    - `f`, a boolean, specifies that matrix entries are symbolic functions
      of the coordinates or are symbolic constants, according to whether
      `f` is True or False.
    - `mode` is a one-letter string.  When`mode` is 'g', 's', or 'a' the
      matrix will be general, symmetric, or antisymmetric.
    """

    def general_matrix(kernel, coords=None, f=False):
        """Returns a general square matrix.  The row index of each entry
        appears as a superscript, while the column index appears as a
        subscript."""
        n = len(coords)
        # Create matrix entries and store in appropriate locations in `G`:
        G = zeros(n, n)
        if f:        # entries are symbolic functions
            for i in range(n):
                for j in range(n):
                    entry = '{' + kernel + '__' + str(coords[i]) + '}_' + str(coords[j])
                    G[i, j] = Function(entry)(*coords)
        else:        # entries are symbolic constants
            for i in range(n):
                for j in range(n):
                    entry = '{' + kernel + '__' + str(coords[i]) + '}_' + str(coords[j])
                    G[i, j] = Symbol(entry, real=True)
        return G

    def symmetric_matrix(kernel, coords=None, f=False):
        """Returns a symmetric matrix.  Entries have a single index, which
        appears as a subscript."""
        n = len(coords)
        # Create and temporarily store matrix entries in `parameters`
        parameters = []
        if f:        # entries are symbolic functions
            for i in range((n*(n+1)//2), 0, -1):
                parameters.append(Function(kernel + '_' + str(i))(*coords))
        else:        # entries are symbolic constants
            for i in range((n*(n+1)//2), 0, -1):
                parameters.append(Symbol(kernel + '_' + str(i), real=True))
        # Transfer entries to symmetric matrix `S`.
        S = zeros(n, n)
        for i in range(n):
            for j in range(i, n):
                S[i, j] = parameters.pop()
                S[j, i] = S[i, j]
        return S

    def antisymmetric_matrix(kernel, coords=None, f=False):
        """Returns an antisymmetric matrix.  Entries have a a single index,
        which appears as a subscript."""
        n = len(coords)
        # Create and temporarily store matrix entries in `parameters`
        parameters = []
        if f:        # entries are symbolic functions
            for i in range((n*(n-1)//2), 0, -1):
                parameters.append(Function(kernel + '_' + str(i))(*coords))
        else:        # entries are symbolic constants
            for i in range((n*(n-1)//2), 0, -1):  # each parameter is a symbol
                parameters.append(Symbol(kernel + '_' + str(i), real=True))
        # Transfer entries to antisymmetric matrix `A`.
        A = zeros(n, n)
        for i in range(n):
            for j in range(i+1, n):
                A[i, j] = parameters.pop()
                A[j, i] = - A[i, j]
        return A

    # Check legitimacy of parameter values:
    if not isinstance(coords, (list, tuple)):
        raise ValueError('coords = ' + str(coords) + ' in Symbolic_Matrix')
    if mode not in ['g', 's', 'a']:
        raise ValueError('mode = ' + str(mode) + ' in Symbolic_Matrix')
    if mode == 'g':
        return general_matrix(kernel, coords, f)
    if mode == 's':
        return symmetric_matrix(kernel, coords, f)
    if mode == 'a':
        return antisymmetric_matrix(kernel, coords, f)
# ## GSG code ends ###


def Matrix_to_dictionary(mat_rep, basis):
    """ Convert matrix representation of linear transformation to dictionary """
    n = len(basis)
    if mat_rep.rows != n or mat_rep.cols != n:
        raise ValueError('Matrix and Basis dimensions not equal for Matrix = ' + str(mat_rep))
    n_range = list(range(n))
    return {
        basis[col]: sum(
            (mat_rep[row, col]*basis[row] for row in n_range), S.Zero
        )
        for col in n_range
    }


# ## GSG code starts ###
def Dictionary_to_Matrix(dict_rep, ga):
    """Returns the matrix representation of that linear transformation on
    geometric algebra ga which has dictionary representation dict_rep."""
    # columns[j] is a list of the entries in the matrix's jth column.
    # columns[j][i] is the (i,j)th entry in the matrix.
    # Matrix[columns] instantiates the transpose of the desired matrix.
    columns = []
    for b in ga.basis:             # b is a basis symbol for ga.
        column = ga.n * [S.Zero]   # Initialize column for dict_rep value at b.
        dict_value = dict_rep[b]   # dict_rep's value at b
        if isinstance(dict_value, mv.Mv):
            dict_value = dict_value.obj
        if dict_value is not S.Zero:
            for coef, base in metric.linear_expand_terms(dict_value):
                row_index = ga.basis.index(base)
                column[row_index] = coef
        columns.append(column)
    return Transpose(Matrix(columns)).doit()
# ## GSG code ends ###


class Lt(printer.GaPrintable):
    r"""
    A Linear Transformation

    Except for the versor representation, the linear transformation
    is stored as a dictionary with basis vector keys and vector
    values ``self.lt_dict`` so that a is a vector :math:`a = a^{i}e_{i}` then

    .. math::
        \mathtt{self(}a\mathtt{)}
            = a^{i} * \mathtt{self.lt\_dict[}e_{i}\mathtt{]}.

    For the versor representation, the linear transformation is
    stored as a versor ``self.V`` so that if a is a
    vector::

        self(a) = self.V.g_invol() * a * self.V.inv()

    where ``self.V.g_invol()`` is the grade involute of ``self.V``.

    Attributes
    ----------
    lt_dict : dict
        the keys are the basis symbols, :math:`e_i`, and the dictionary
        entries are the object vector images (linear combination of sympy
        non-commutative basis symbols) of the keys so that if ``L`` is the
        linear transformation then::

            L(e_i) = self.Ga.mv(L.lt_dict[e_i])

    """

    @staticmethod
    def setup(ga):
        # galgebra 0.5.0
        warnings.warn(
            "Lt.setup(ga) is deprecated, use `ga.coords` and `ga.coord_vec` "
            "directly.", DeprecationWarning, stacklevel=2)
        return ga.coords, ga.coord_vec

    @property
    def coords(self):
        # galgebra 0.6.0
        warnings.warn(
            "lt.coords is deprecated, use `lt.Ga.coords` instead.",
            DeprecationWarning, stacklevel=2)
        return self.Ga.coords

    @property
    def X(self):
        # galgebra 0.6.0
        warnings.warn(
            "lt.X is deprecated, use `lt.Ga.coord_vec` instead.",
            DeprecationWarning, stacklevel=2)
        return self.Ga.coord_vec

    @property
    def mode(self):
        # galgebra 0.6.0
        warnings.warn(
            "lt.mode is deprecated, inspect lt.matrix() and its transpose to "
            "determine symmetry",
            DeprecationWarning, stacklevel=2)
        m = self.matrix()
        if m == m.T:
            return 's'
        elif m == -m.T:
            return 'a'
        else:
            return 'g'

    @property
    def fct_flg(self):
        # galgebra 0.6.0
        warnings.warn(
            "lt.fct_flg is deprecated, inspect lt.matrix().free_symbols to "
            "determine coordinate-dependence",
            DeprecationWarning, stacklevel=2)
        if self.Ga.coords is None:
            return False
        return set(self.Ga.coords) <= self.matrix().free_symbols

    def __init__(self, *args, ga, f=False, mode='g'):
        """
        __init__(self, *args, ga, **kwargs)

        Note this constructor is overloaded, based on the type of the
        positional argument:

        .. class:: Lt(lt_dict: Dict[Expr, Expr], /, *, ga)
            :noindex:

            Construct from a dictionary mapping source basis blade expressions
            to multivectors.

        .. class:: Lt(lt_matrix: Matrix, /, *, ga)
            :noindex:

            Construct from the operation of matrix pre-multiplication.

        # ## GSG code starts ###
        .. class:: Lt(lt_list: list, /, *, ga)
            :noindex:

            Construct from a list of lists, the j_th list of which contains
            the coefficients of j_th image vector's basis expansion.
        # ## GSG code ends ###

        .. class:: Lt(versor: mv.Mv, /, *, ga)
            :noindex:

            Construct from a not-necessarily-normalized versor.

        .. class:: Lt(func: Callable[[mv.Mv], mv.Mv], /, *, ga)
            :noindex:

            Construct from a function, which is tested for linearity.

        .. class:: Lt(s: str, /, *, ga, f=False, mode='g')
            :noindex:

            Construct an appropriate matrix from a string `s`.


        Parameters
        ----------
        ga : Ga
            Geometric algebra which is both domain and codomain of this transformation
        f : bool
            True if Lt if function of coordinates. Only supported in the string
            constructor
        mode : str
            g:general, s:symmetric, a:antisymmetric transformation.
            Only supported in the string constructor.
        """

        mat_rep = args[0]
        self.Ga = ga
        self.lt_dict = {}
        self.mat = None
        self.versor = False
        # self.V, self.Vrev, and self.Vqform are never actually used in the current
        # implementation of orthogonal outermorphisms created by a versor input.
        self.V = None
        self.Vrev = None
        self.Vqform = None

        if isinstance(mat_rep, dict):         # Dictionary input
            for key in mat_rep:
                self.lt_dict[key] = mat_rep[key]

        elif isinstance(mat_rep, list):       # List input
            if not isinstance(mat_rep[0], list):
                # At this point mat_rep[i] is the desired image vector for the
                # i_th basis image vectors.
                for lt_i, base in zip(mat_rep, self.Ga.basis):
                    self.lt_dict[base] = sym(lt_i)
            else:
                # mat_rep = map(list, zip(*mat_rep))  # Transpose list of lists
                for row, base1 in zip(mat_rep, self.Ga.basis):
                    tmp = 0
                    for col, base2 in zip(row, self.Ga.basis):
                        tmp += col * base2
                    self.lt_dict[base1] = tmp

        # ## GSG code starts ###
        elif isinstance(mat_rep, Matrix):    # Matrix input
            self.lt_dict = Matrix_to_dictionary(mat_rep, self.Ga.basis)
        # ## GSG code ends ###

        # ## GSG code starts. ###
        # This code segment uses versor `mat_rep` and a sandwich product only to
        # create a linear vector-valued function of vector. That function is then
        # used to create a dictionary-based outermorphism. Evaluation of the
        # outermorphism on a multivector is by dictionary lookup, not by a
        # sandwich product of the multivector with the versor.
        elif isinstance(mat_rep, mv.Mv):     # Versor input
            if not mat_rep.is_versor:
                raise ValueError(mat_rep, 'is not a versor in Versor input for Lt!\n')
            V = mat_rep
            Vg_invol = V.g_invol()
            Vinv = V.inv()
            outermorphism = ga.lt(lambda x: Vg_invol * x * Vinv)
            self.lt_dict = simplify(outermorphism.lt_dict)
        # ## GSG code ends ###

        # ## GSG code starts ###
        elif isinstance(mat_rep, str):    # (One-letter) string input
            Amat = Symbolic_Matrix(mat_rep, coords=self.Ga.coords, f=f, mode=mode)
            if mode == 'g':
                self.__init__(Amat, ga=self.Ga)
            elif mode in ['s', 'a']:
                self.__init__(self.Ga.g_inv * Amat, ga=self.Ga)
        # ## GSG code ends ###

        elif callable(mat_rep):    # Linear multivector function input
            # Function is tested for linearity before use.
            F = mat_rep
            a = mv.Mv('a', 'vector', ga=self.Ga)
            b = mv.Mv('b', 'vector', ga=self.Ga)
            if F(a + b) != F(a) + F(b):
                raise ValueError('{} is not linear'.format(F))
            self.lt_dict = {}
            for base in self.Ga.basis:
                out = F(mv.Mv(base, ga=self.Ga))
                if not out.is_vector():
                    raise ValueError('{} must return vectors'.format(F))
                self.lt_dict[base] = out.obj
        else:
            raise TypeError("Unsupported argument type {}".format(type(mat_rep)))

    @_cached_property
    def mv_dict(self) -> Mapping[Expr, Expr]:
        # dict for linear transformation of multivector
        if self.versor:
            # no lt_dict
            return None
        return {
            blade: reduce(
                self.Ga.wedge,
                (self.Ga.basis[i].xreplace(self.lt_dict) for i in index),
                S.One
            )
            for index, blade in self.Ga.indexes_to_blades_dict.items()
        }

    def __call__(self, M, obj=False):
        r"""
        Returns the image of multivector :math:`M` under the linear transformation
        :math:`L`. :math:`{{L}\lp{M}\rp}` is defined by
        the linearity of :math:`L`,
        the vector values :math:`{{L}\lp{{{\eb}}_{j}}\rp }`, and the definition
        :math:`{{L}\lp{{{\eb}}_{j_{1}}{\wedge}\dots{\wedge}{{\eb}}_{j_{r}}}\rp}={{L}\lp{{{\eb}}_{j_{1}}}\rp}{\wedge}\dots{\wedge}{{L}\lp{{{\eb}}_{j_{r}}}\rp}`.
        """

        if isinstance(M, mv.Mv) and self.Ga != M.Ga:
            raise ValueError('In Lt call Lt and argument refer to different vector spaces')

        # ## GSG code starts ###
        # Given the current way an outermorphism is created from a versor input,
        # self.versor will always be false; hence the following code fragment will
        # never execute.
        if self.versor:
            # Sandwich M or M's grade involute depending on whether versor self.V
            # is even or odd.
            if self.V == self.V.odd():
                V_M_Vrev = self.V * M.g_invol() * self.Vrev
            elif self.V == self.V.even():
                V_M_Vrev = self.V * M * self.Vrev
            else:
                raise ValueError('self.V is not a versor in  __call__')
            # Divide by normalization factor self.Vqform to convert sandwiching
            # between self.V and its reverse to sandwiching between self.V and
            # its inverse.
            V_M_Vinv = 1/(self.Vqform) * V_M_Vrev
            if obj:
                return V_M_Vinv.obj
            else:
                return V_M_Vinv
        # ## GSG code ends ###

        if isinstance(M, mv.Mv):
            if M.is_vector():
                lt_M = M.obj.xreplace(self.lt_dict)
                if obj:
                    return lt_M
                else:
                    return mv.Mv(lt_M, ga=self.Ga)
            else:
                mv_obj = M.obj
        else:
            mv_obj = mv.Mv(M, ga=self.Ga).obj

        lt_M = mv_obj.xreplace(self.mv_dict)
        if obj:
            return lt_M
        else:
            return mv.Mv(lt_M, ga=self.Ga)

    def __add__(self, LT):

        if self.Ga != LT.Ga:
            raise ValueError("Attempting addition of Lt's from different geometric algebras")

        self_add_LT = copy(self.lt_dict)
        for key in list(LT.lt_dict.keys()):
            if key in self_add_LT:
                self_add_LT[key] = metric.collect(self_add_LT[key] + LT.lt_dict[key], self.Ga.basis)
            else:
                self_add_LT[key] = LT.lt_dict[key]
        return Lt(self_add_LT, ga=self.Ga)

    def __sub__(self, LT):

        if self.Ga != LT.Ga:
            raise ValueError("Attempting subtraction of Lt's from different geometric algebras")

        self_add_LT = copy(self.lt_dict)
        for key in list(LT.lt_dict.keys()):
            if key in self_add_LT:
                self_add_LT[key] = metric.collect(self_add_LT[key] - LT.lt_dict[key], self.Ga.basis)
            else:
                self_add_LT[key] = -LT.lt_dict[key]
        return Lt(self_add_LT, ga=self.Ga)

    def __mul__(self, LT):

        if isinstance(LT, Lt):

            if self.Ga != LT.Ga:
                raise ValueError("Attempting multiplication of Lt's from different geometric algebras")
            self_mul_LT = {}
            for base in LT.lt_dict:
                self_mul_LT[base] = self(LT(base, obj=True), obj=True)
            for key in self_mul_LT:
                self_mul_LT[key] = metric.collect(expand(self_mul_LT[key]), self.Ga.basis)
            return Lt(self_mul_LT, ga=self.Ga)
        else:
            self_mul_LT = {}
            for key in self.lt_dict:
                self_mul_LT[key] = LT * self.lt_dict[key]
            return Lt(self_mul_LT, ga=self.Ga)

    def __rmul__(self, LT):

        if not isinstance(LT, Lt):
            self_mul_LT = {}
            for key in self.lt_dict:
                self_mul_LT[key] = LT * self.lt_dict[key]
            return Lt(self_mul_LT, ga=self.Ga)
        else:
            raise TypeError('Cannot have LT as left argument in Lt __rmul__\n')

    # ## GSG code starts ###
    def det(self) -> Expr:    # det(L) defined by L(E) = det(L)E
        r"""
        - Returns the determinant of the linear transformation :math:`L`,
          defined by :math:`\det(L) = L(E) E^{-1}`, where :math:`E` is the
          basis blade for the pseudoscalar grade space.
        - Expression returned is a real SymPy scalar, not a GAlgebra 0-vector.
        """
        return (self(self.Ga.e) * self.Ga.e.inv()).scalar()
    # ## GSG code ends ###

    def tr(self) -> Expr:  # tr(L) defined by tr(L) = grad|L(x)
        r"""
        Returns the trace (a scalar) of the linear transformation,
        :math:`L`, defined by :math:`{{\operatorname{tr}}\lp {L} \rp }=\nabla_{a}\cdot{{L}\lp {a} \rp }`
        where :math:`a` is a vector in the tangent space.
        """
        connect_flg = self.Ga.connect_flg
        self.Ga.connect_flg = False
        F_x = mv.Mv(self(self.Ga.coord_vec, obj=True), ga=self.Ga)
        tr_F = (self.Ga.grad | F_x).scalar()
        self.Ga.connect_flg = connect_flg
        return tr_F

    r'''
    def adj(self) -> 'Lt':
        r"""
        Returns the adjoint :math:`{\bar{L}}`(a linear transformation) of linear
        transformation :math:`L`, defined by
        :math:`a\cdot{{L}\lp {b} \rp } = b\cdot{{\bar{L}}\lp {a} \rp }`
        where :math:`a` and :math:`b` are any two vectors in the tangent space.
        """
        self_adj = []
        for e_j in self.Ga.basis:
            s = S.Zero
            for e_i, er_i in zip(self.Ga.basis, self.Ga.r_basis):
                s += er_i * self.Ga.hestenes_dot(e_j, self(e_i, obj=True))
            if self.Ga.is_ortho:
                self_adj.append(expand(s))
            else:
                self_adj.append(expand(s) / self.Ga.e_sq)
        return Lt(self_adj, ga=self.Ga)
    '''

    # ## GSG code starts ###
    def adj(self) -> 'Lt':
        r"""
        Returns the adjoint transformation :math:`{\bar{L}}` of linear
        transformation :math:`L`, defined by
        :math:`a\cdot{{L}\lp {b} \rp } = b\cdot{{\bar{L}}\lp {a} \rp }`,
        where :math:`a` and :math:`b` are any two vectors in the tangent space.
        """
        matrix_of_adjoint = self.Ga.g_inv * self.matrix().T * self.Ga.g
        return self.Ga.lt(matrix_of_adjoint)
    # ## GSG code ends ###

    # ## GSG code starts ###
    def is_singular(self):
        """Returns `True` if and only if  linear transformation `self` is singular."""
        E = self.Ga.E()
        return simplify((self(E) < E.inv()).scalar()) == S.Zero
    # ## GSG code ends

    # ## GSG code starts ###
    def inv(self):
        """Returns compositional inverse of linear transformation`self`.
        Assumes transformation is nonsingular.  If `self` is a versor based
        transformation, its inverse will also be versor based."""
        if self.versor:
            return self.Ga.lt(self.V.rev())
        if not self.is_singular():
            return self.Ga.lt(Matrix(self.matrix().inv()))
        else:
            raise ValueError('transformation in inv() is non-invertible')
    # ## GSG code ends ###

    def _sympystr(self, print_obj):

        if self.versor:    # ## GSG: changed `self.spinor` to `self.versor` ###
            return 'R = ' + print_obj._print(self.V)
        else:
            pre = 'Lt('
            s = ''
            for base in self.Ga.basis:
                if base in self.lt_dict:
                    s += pre + print_obj._print(base) + ') = ' + print_obj._print(mv.Mv(self.lt_dict[base], ga=self.Ga)) + '\n'
                else:
                    s += pre + print_obj._print(base) + ') = 0\n'
            return s[:-1]

    # ## GSG code starts ###
    def _latex(self, print_obj):
        parts = []
        for base in self.Ga.basis:           # base is a basis symbol
            if self.versor:
                b = mv.Mv(base, ga=self.Ga)  # b is the corresponding basis vector
                if self.V == self.V.odd():
                    unnormalized_image = self.V * (b.g_invol()) * self.Vrev
                elif self.V == self.V.even():
                    unnormalized_image = self.V * b * self.Vrev
                else:
                    raise ValueError('self.V is not a versor in  _latex')
                image = 1/(self.Vqform) * unnormalized_image
            else:
                image = mv.Mv(self.lt_dict.get(base, S.Zero), ga=self.Ga)
            parts.append(print_obj._print(base) + ' &\\mapsto ' + print_obj._print(image))
        return '\\left\\{ \\begin{aligned} ' + ' \\\\ '.join(parts) + ' \\end{aligned} \\right\\}'
    # ## GSG code ends ###

    def Fmt(self, fmt=1, title=None) -> printer.GaPrintable:
        return printer._FmtResult(self, title)

    # ## GSG code starts ###
    def matrix(self) -> Matrix:
        r"""
        Returns the matrix :math:`[{L__i}_j]` defined for linear transformation
        :math:`L` by :math:`L({\eb}_j)=\sum_i {L__i}_j \eb}_i`.
        """
        if self.mat is not None:
            return self.mat
        elif self.versor:
            self.lt_dict = {}
            for base in self.Ga.basis:
                self.lt_dict[base] = self(base).simplify()
            self.versor = False    # temporary change of self.versor
            mat = self.matrix()
            self.versor = True     # reverse change to self.versor
            return mat
        else:
            self.mat = Dictionary_to_Matrix(self.lt_dict, self.Ga)
            return self.mat.doit()
    # ## GSG code ends ###


class Mlt(printer.GaPrintable):
    r"""
    A multilinear transformation (mlt) is a multilinear multivector function of
    a list of vectors (``*args``) :math:`F(v_1,...,v_r)` where for any argument slot
    :math:`j` we have (:math:`a` is a scalar and :math:`u_j` a vector)

    .. math::
          F(v_1,...,a*v_j,...,v_r) &= a*F(v_1,...,v_j,...,v_r) \\
          F(v_1,...,v_j+u_j,...,v_r) &= F(v_1,...,v_j,...,v_r) + F(v_1,...,u_j,...,v_r).

    If F and G are two :class:`Mlt`\ s with the same number of argument slots then the sum is

    .. math:: (F+G)F(v_1,...,v_r) = F(v_1,...,v_r) + G(v_1,...,v_r).

    If :math:`F` and :math:`G` are two :class:`Mlt`\ s with :math:`r` and :math:`s`
    argument slots then their product is

    .. math:: (F*G)(v_1,...,v_r,...,v_{r+s}) = F(v_1,...,v_r)*G(v_{r+1},...,v_{r+s}),

    where :math:`*` is any of the multivector multiplicative operations.
    The derivative of a :class:`Mlt` with is defined as the directional derivative with respect
    to the coordinate vector (we assume :math:`F` is implicitely a function of the
    coordinates)

    .. math:: F(v_1,...,v_r;v_{r+1}) = (v_{r+1} \bullet \nabla)F(v_1,...,v_j,...,v_r).

    The contraction of a :class:`Mlt` between slots :math:`j` and :math:`k` is defined as the
    geometric derivative of :math:`F` with respect to slot :math:`k` and the inner geometric
    derivative with respect to slot :math:`j` (this gives the standard tensor
    definition of contraction for the case that :math:`F` is a scalar function)

    .. math::

        \operatorname{Contract}(i,j,F)
            &= \nabla_i \bullet (\nabla_j F(v_1,...,v_i,...,v_j,...,v_r)) \\
            &= \nabla_j \bullet (\nabla_i F(v_1,...,v_i,...,v_j,...,v_r)).

    This returns a :class:`Mlt`\ with slot :math:`i` and :math:`j` removed.
    """

    @staticmethod
    def subs(Ga, anew):
        #  Generate coefficient substitution list for new Mlt slot
        #  vectors (arguments) where anew is a list of slot vectors
        #  to be substituted for the old slot vectors.
        #  This is used when one wishes to substitute specific vector
        #  values into the Mlt such as the basis/reciprocal basis vectors.
        sub_lst = []
        for i, a in enumerate(anew):
            acoefs = a.get_coefs(1)
            sub_lst += list(zip(Ga._mlt_pdiffs[i], acoefs))
        return sub_lst

    @staticmethod
    def increment_slots(nargs, Ga):
        # Increment cache of available slots (vector variables) if needed for Mlt class
        n_a = len(Ga._mlt_a)
        if n_a < nargs:
            for i in range(n_a, nargs):
                #  New slot variable with coefficients a_{n_a}__k
                a = Ga.mv('a_' + str(i + 1), 'vector')
                #  Append new slot variable a_j
                Ga._mlt_a.append(a)
                #  Append slot variable coefficients a_j__k for purpose
                #  of differentiation
                coefs = a.get_coefs(1)
                Ga._mlt_pdiffs.append(coefs)
                Ga._mlt_acoefs += coefs

    @staticmethod
    def extact_basis_indexes(Ga):
        # galgebra 0.5.0
        warnings.warn(
            "`Mlt.extact_basis_indexes(ga)` is deprecated, use `ga.basis_super_scripts`",
            DeprecationWarning, stacklevel=2)
        return Ga.basis_super_scripts

    def _sympystr(self, print_obj):
        return print_obj._print(self.fvalue)

    def _latex(self, print_obj):
        if self.nargs <= 1:
            return print_obj._print(self.fvalue)
        expr_lst = Mlt.expand_expr(self.fvalue, self.Ga)
        latex_str = '\\begin{aligned} '
        first = True
        lcnt = print_obj._settings['galgebra_mlt_lcnt']
        cnt = 1  # Component count on line
        for term in expr_lst:
            coef_str = str(term[0])
            coef_latex = print_obj._print(term[0])
            term_add_flg = isinstance(term[0], Add)
            if term_add_flg:
                coef_latex = r'\left ( ' + coef_latex + r'\right ) '
            if first:
                first = False
            else:
                if coef_str[0].strip() != '-' or term_add_flg:
                    coef_latex = ' + ' + coef_latex
            for aij in term[1]:
                coef_latex += print_obj._print(aij) + ' '
            if cnt == 1:
                latex_str += ' & ' + coef_latex
            else:
                latex_str += coef_latex
            if cnt % lcnt == 0:
                latex_str += '\\\\ '
                cnt = 1
            else:
                cnt += 1
        if lcnt == len(expr_lst) or lcnt == 1:
            latex_str = latex_str[:-3]
        latex_str = latex_str + ' \\end{aligned} '
        return latex_str

    def Fmt(self, lcnt=1, title=None) -> printer.GaPrintable:
        """
        Set format for printing of Tensors

        Parameters
        ----------
        lcnt :
            Number of components per line

        Notes
        -----
        Usage for tensor T example is::

            T.fmt('2', 'T')

        output is::

            print 'T = '+str(A)

        with two components per line.  Works for both standard printing and
        for latex.
        """
        obj = printer._WithSettings(self, dict(galgebra_mlt_lcnt=lcnt))
        return printer._FmtResult(obj, title)

    @staticmethod
    def expand_expr(expr, ga):
        lst_expr = []
        expr = expand(expr)
        for term in expr.args:
            coef = S.One
            a_lst = []
            for factor in term.args:
                if factor in ga._mlt_acoefs:
                    a_lst.append(factor)
                else:
                    coef *= factor
            a_lst = tuple([x for x in a_lst if x in ga._mlt_acoefs])
            b_lst = tuple([ga._mlt_acoefs.index(x) for x in a_lst])
            lst_expr.append((coef, a_lst, b_lst))
        lst_expr = sorted(lst_expr, key=lambda x: x[2])
        new_lst_expr = []
        previous = (-1,)
        first = True
        a = None
        for term in lst_expr:
            if previous == term[2]:
                coef += term[0]
                previous = term[2]
            else:
                if not first:
                    new_lst_expr.append((coef, a))
                else:
                    first = False
                coef = term[0]
                previous = term[2]
                a = term[1]
        new_lst_expr.append((coef, a))
        return new_lst_expr

    def __init__(self, f, Ga, nargs=None, fct=False):
        #  f is a function, a multivector, a string, or a component expression
        #  self.f is a function or None such as T | a_1 where T and a_1 are vectors
        #  self.fvalue is a component expression such as
        #  T_x*a_1__x+T_y*a_1__y+T_z*a_1__z for a rank 1 tensor in 3 space and all
        #  symbols are sympy real scalar symbols
        self.Ga = Ga
        if isinstance(f, mv.Mv):
            if f.is_vector():  # f is vector T = f | a1
                self.nargs = 1
                Mlt.increment_slots(1, Ga)
                self.fvalue = (f | Ga._mlt_a[0]).obj
                self.f = None
            else:  # To be inplemented for f a general pure grade mulitvector
                self.nargs = nargs
                self.fvalue = f
                self.f = None
        elif isinstance(f, Lt):  # f is linear transformation T = a1 | f(a2)
            self.nargs = 2
            Mlt.increment_slots(2, Ga)
            self.fvalue = (Ga._mlt_a[0] | f(Ga._mlt_a[1])).obj
            self.f = None
        elif isinstance(f, str) and nargs is not None:
            self.f = None
            self.nargs = nargs
            Mlt.increment_slots(nargs, Ga)
            self.fvalue = S.Zero
            for t_index, a_prod in zip(itertools.product(self.Ga.basis_super_scripts, repeat=self.nargs),
                                       itertools.product(*self.Ga._mlt_pdiffs)):
                name = '{}_{}'.format(f, ''.join(map(str, t_index)))
                if fct:  # Tensor field
                    coef = Function(name, real=True)(*self.Ga.coords)
                else:  # Constant Tensor
                    coef = symbols(name, real=True)
                self.fvalue += reduce(lambda x, y: x*y, a_prod, coef)

        else:
            if isinstance(f, types.FunctionType):  # Tensor defined by general multi-linear function
                args = inspect.getfullargspec(f)[0]
                self.nargs = len(args)
                self.f = f
                Mlt.increment_slots(self.nargs, Ga)
                self.fvalue = f(*tuple(Ga._mlt_a[0:self.nargs]))
            else:  # Tensor defined by component expression
                raise NotImplementedError
                # self.f = None
                # self.nargs = len(args)  # args isn't defined, which is why we raise NotImplementedError
                # Mlt.increment_slots(self.nargs, Ga)
                # self.fvalue = f

    def __call__(self, *args):
        """
        Evaluate the multilinear function for the given vector arguments.
        Note that a sympy scalar is returned, *not* a multilinear function.
        """
        if len(args) == 0:
            return self.fvalue
        if self.f is not None:
            return self.f(*args)
        else:
            sub_lst = []
            for x, ai in zip(args, self.Ga._mlt_pdiffs):
                for r_base, aij in zip(self.Ga.r_basis_mv, ai):
                    sub_lst.append((aij, (r_base | x).scalar()))
            return self.fvalue.subs(sub_lst, simultaneous=True)

    def __add__(self, X):
        if isinstance(X, Mlt):
            if self.nargs == X.nargs:
                return Mlt(self.fvalue + X.fvalue, self.Ga, self.nargs)
            else:
                raise ValueError('In Mlt add number of args not the same\n')
        else:
            raise TypeError('In Mlt add second argument not an Mkt\n')

    def __sub__(self, X):
        if isinstance(X, Mlt):
            if self.nargs == X.nargs:
                return Mlt(self.fvalue - X.fvalue, self.Ga, self.nargs)
            else:
                raise ValueError('In Mlt sub number of args not the same\n')
        else:
            raise TypeError('In Mlt sub second argument not an Mlt\n')

    def __mul__(self, X):
        if isinstance(X, Mlt):
            nargs = self.nargs + X.nargs
            Mlt.increment_slots(nargs, self.Ga)
            self_args = self.Ga._mlt_a[:self.nargs]
            X_args = X.Ga._mlt_a[self.nargs:nargs]
            value = (self(*self_args) * X(*X_args)).expand()
            return Mlt(value, self.Ga, nargs)
        else:
            return Mlt(X * self.fvalue, self.Ga, self.nargs)

    def __xor__(self, X):
        if isinstance(X, Mlt):
            nargs = self.nargs + X.nargs
            Mlt.increment_slots(nargs, self.Ga)
            value = self(*self.Ga._mlt_a[:self.nargs]) ^ X(*X.Ga._mlt_a[self.nargs:nargs])
            return Mlt(value, self.Ga, nargs)
        else:
            return Mlt(X * self.fvalue, self.Ga, self.nargs)

    def __or__(self, X):
        if isinstance(X, Mlt):
            nargs = self.nargs + X.nargs
            Mlt.increment_slots(nargs, self.Ga)
            value = self(*self.Ga._mlt_a[:self.nargs]) | X(*X.Ga._mlt_a[self.nargs:nargs])
            return Mlt(value, self.Ga, nargs)
        else:
            return Mlt(X * self.fvalue, self.Ga, self.nargs)

    def dd(self):
        Mlt.increment_slots(self.nargs + 1, self.Ga)
        dd_fvalue = (self.Ga._mlt_a[self.nargs] | self.Ga.grad) * self.fvalue
        return Mlt(dd_fvalue, self.Ga, self.nargs + 1)

    def pdiff(self, slot: int):
        r"""
        Returns gradient of tensor, ``T``, with respect to slot vector.

        For example if the tensor is :math:`{{T}\lp {a_{1},a_{2}} \rp }` then ``T.pdiff(2)`` is :math:`\nabla_{a_{2}}T`. Since ``T`` is a scalar function,
        ``T.pdiff(2)`` is a vector function.
        """
        # Take geometric derivative of mlt with respect to slot argument
        self.Ga.dslot = slot - 1
        return self.Ga.grad * self.Ga.mv(self.fvalue)

    @staticmethod
    def remove_slot(mv, slot, nargs, ga):
        if slot == nargs:
            return mv
        for islot in range(slot, nargs):
            mv = mv.subs(list(zip(ga._mlt_pdiffs[islot], ga._mlt_pdiffs[islot - 1])))
        return mv

    def contract(self, slot1: int, slot2: int):
        """
        Returns contraction of tensor between ``slot1`` and ``slot2`` where
        ``slot1`` is the index of the first vector argument and ``slot2`` is the
        index of the second vector argument of the tensor.

        For example if we have a rank two tensor, ``T(a1, a2)``, then
        ``T.contract(1, 2)`` is the contraction of ``T``.
        For this case since there are only two slots, there can only be one
        contraction.
        """
        min_slot = min(slot1, slot2)
        max_slot = max(slot1, slot2)
        cnargs = self.nargs - 2
        self.Ga.dslot = min_slot - 1
        grad_self = self.Ga.grad * self.Ga.mv(self.fvalue)
        grad_self = Mlt.remove_slot(grad_self.obj, min_slot, self.nargs, self.Ga)
        self.Ga.dslot = max_slot - 2
        div_grad_self = self.Ga.grad | self.Ga.mv(grad_self)
        div_grad_self = Mlt.remove_slot(div_grad_self.obj, max_slot - 1, self.nargs - 1, self.Ga)
        return Mlt(div_grad_self, self.Ga, cnargs)

    def cderiv(self):
        """
        Returns covariant derivative of tensor field.

        If ``T`` is a tensor of rank :math:`k` then ``T.cderiv()`` is a tensor
        of rank :math:`k+1`. The operation performed is defined in section
        :ref:`MLtrans`.
        """
        Mlt.increment_slots(self.nargs + 1, self.Ga)
        agrad = self.Ga._mlt_a[self.nargs] | self.Ga.grad
        CD = Mlt((agrad * self.Ga.mv(self.fvalue)).obj, self.Ga, self.nargs + 1)
        if CD != 0:
            CD = CD.fvalue
        for i in range(self.nargs):
            args = self.Ga._mlt_a[:self.nargs]
            tmp = agrad * self.Ga._mlt_a[i]
            if tmp.obj != 0:
                args[i] = tmp
                CD = CD - self(*args)
        CD = Mlt(CD, self.Ga, self.nargs + 1)
        return CD

    def expand(self):
        self.fvalue = expand(self.fvalue)
        return self

    def comps(self):
        basis = self.Ga.mv()
        rank = self.nargs
        ndim = len(basis)
        i_indexes = itertools.product(list(range(ndim)), repeat=rank)
        indexes = itertools.product(basis, repeat=rank)
        output = ''
        for i, (e, i_index) in enumerate(zip(indexes, i_indexes)):
            if i_index[-1] % ndim == 0:
                print('')
            output += str(i)+':'+str(i_index)+':'+str(self(*e)) + '\n'
        return output


# ## GSG code starts ###
def det(L: Lt) -> Expr:    # det(L) defined by L(E) = det(L)E
    r"""
    - Returns the determinant of the linear transformation :math:`L`,
      defined by :math:`\det(L) = L(E) E^{-1}`, where :math:`E` is the
      basis blade for the pseudoscalar grade space.
    - Expression returned is a real SymPy scalar, not a GAlgebra 0-vector.
    """
    return L.det()
# ## GSG code ends ###


# ## GSG code starts ###
def sym(v):
    """
    Returns that linear combination of basis vector symbols which corresponds
    to vector v, itself a linear combination of basis vectors.
    """
    # Obtain the coefficients in basis vector expansion of `v`.
    # Then construct and return corresponding basis vector symbol expansion.
    coefs = v.blade_coefs(v.Ga.mv())
    return sum(coefs[j]*v.Ga.basis[j] for j in range(v.Ga.n))
# ## GSG code ends ###
