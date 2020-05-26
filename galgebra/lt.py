"""
Multivector Linear Transformation
"""

import inspect
import types
import itertools
import warnings
from copy import copy
from functools import reduce

from sympy import (
    expand, symbols, Matrix, Transpose, zeros, Symbol, Function, S, Add, Expr
)
from sympy.printing.latex import LatexPrinter as _LatexPrinter
from sympy.printing.str import StrPrinter as _StrPrinter


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


def Symbolic_Matrix(root, coords=None, mode='g', f=False, sub=True):
    if sub:
        pos = '_'
    else:
        pos = '__'
    if isinstance(coords, (list, tuple)):
        n = len(coords)
        n_range = range(n)
        mat = zeros(n)
        if mode == 'g':  # General symbolic matrix
            for row in n_range:
                row_index = str(coords[row])
                for col in n_range:
                    col_index = str(coords[col])
                    element = root + pos + row_index + col_index
                    if not f:
                        mat[row, col] = Symbol(element, real=True)
                    else:
                        mat[row, col] = Function(element)(*coords)

        elif mode == 's':  # Symmetric symbolic matrix
            for row in n_range:
                row_index = str(coords[row])
                for col in n_range:
                    col_index = str(coords[col])
                    if row <= col:
                        element = root + pos + row_index + col_index
                    else:
                        element = root + pos + col_index + row_index
                    if not f:
                        mat[row, col] = Symbol(element, real=True)
                    else:
                        mat[row, col] = Function(element)(*coords)

        elif mode == 'a':  # Asymmetric symbolic matrix
            for row in n_range:
                row_index = str(coords[row])
                for col in n_range:
                    col_index = str(coords[col])
                    if row <= col:
                        sign = S(1)
                        element = root + pos + row_index + col_index
                    else:
                        sign = -S(1)
                        element = root + pos + col_index + row_index
                    if row == col:
                        sign = S(0)
                    if not f:
                        mat[row, col] = sign * Symbol(element, real=True)
                    else:
                        mat[row, col] = sign * Function(element)(*coords)
        else:
            raise ValueError('In Symbolic_Matrix mode = ' + str(mode))
    else:
        raise ValueError('In Symbolic_Matrix coords = ' + str(coords))
    return mat


def Matrix_to_dictionary(mat_rep, basis):
    """ Convert matrix representation of linear transformation to dictionary """
    dict_rep = {}
    n = len(basis)
    if mat_rep.rows != n or mat_rep.cols != n:
        raise ValueError('Matrix and Basis dimensions not equal for Matrix = ' + str(mat_rep))
    n_range = list(range(n))
    for row in n_range:
        dict_rep[basis[row]] = S(0)
        for col in n_range:
            dict_rep[basis[row]] += mat_rep[col, row]*basis[col]
    return dict_rep


def Dictionary_to_Matrix(dict_rep, ga):
    """ Convert dictionary representation of linear transformation to matrix """
    basis = list(dict_rep.keys())
    n = len(basis)
    n_range = list(range(n))
    lst_mat = []  # list representation of sympy matrix
    for row in n_range:
        e_row = ga.basis[row]
        lst_mat_row = n * [S(0)]

        if e_row in basis:  # If not in basis row all zeros
            element = dict_rep[e_row]
            if isinstance(element, mv.Mv):
                element = element.obj
            coefs, bases = metric.linear_expand(element)
            for coef, base in zip(coefs, bases):
                index = ga.basis.index(base)
                lst_mat_row[index] = coef

        lst_mat.append(lst_mat_row)
    return Transpose(Matrix(lst_mat))


class Lt(printer.GaPrintable):
    r"""
    A Linear Transformation

    Except for the spinor representation the linear transformation
    is stored as a dictionary with basis vector keys and vector
    values ``self.lt_dict`` so that a is a vector :math:`a = a^{i}e_{i}` then

    .. math::
        \mathtt{self(}a\mathtt{)}
            = a^{i} * \mathtt{self.lt\_dict[}e_{i}\mathtt{]}.

    For the spinor representation the linear transformation is
    stored as the even multivector ``self.R`` so that if a is a
    vector::

        self(a) = self.R * a * self.R.rev().

    Attributes
    ----------
    lt_dict : dict
        the keys are the basis symbols, :math:`e_i`, and the dictionary
        entries are the object vector images (linear combination of sympy
        non-commutative basis symbols) of the keys so that if ``L`` is the
        linear transformation then::

            L(e_i) = self.Ga.mv(L.lt_dict[e_i])

    """

    mat_fmt = False

    @staticmethod
    def setup(ga):
        return ga.coords, ga.coord_vec

    @staticmethod
    def format(mat_fmt=False):
        Lt.mat_fmt = mat_fmt

    def __init__(self, *args, ga, f=False, mode='g'):
        """
        Parameters
        ----------
        ga :
            Name of metric (geometric algebra)
        f : bool
            True if Lt if function of coordinates
        mode : str
            g:general, s:symmetric, a:antisymmetric transformation
        """
        mat_rep = args[0]
        self.fct_flg = f
        self.mode = mode
        self.Ga = ga
        self.coords = ga.coords
        self.X = ga.coord_vec
        self.spinor = False
        self.rho_sq = None

        self.lt_dict = {}
        self.mv_dict = None
        self.mat = None

        if isinstance(mat_rep, tuple):  # tuple input
            for key in mat_rep:
                self.lt_dict[key] = mat_rep[key]

        elif isinstance(mat_rep, dict):  # Dictionary input
            for key in mat_rep:
                self.lt_dict[key] = mat_rep[key]

        elif isinstance(mat_rep, list):  # List of lists input
            if not isinstance(mat_rep[0], list):
                for lt_i, base in zip(mat_rep, self.Ga.basis):
                    self.lt_dict[base] = lt_i
            else:
                # mat_rep = map(list, zip(*mat_rep))  # Transpose list of lists
                for row, base1 in zip(mat_rep, self.Ga.basis):
                    tmp = 0
                    for col, base2 in zip(row, self.Ga.basis):
                        tmp += col * base2
                    self.lt_dict[base1] = tmp

        elif isinstance(mat_rep, Matrix):  # Matrix input
            self.mat = mat_rep
            mat_rep = self.mat * self.Ga.g_inv
            self.lt_dict = Matrix_to_dictionary(mat_rep, self.Ga.basis)

        elif isinstance(mat_rep, mv.Mv):  # Spinor input
            self.spinor = True
            self.R = mat_rep
            self.Rrev = mat_rep.rev()
            self.rho_sq = self.R * self.Rrev
            if self.rho_sq.is_scalar():
                self.rho_sq = self.rho_sq.scalar()
                if self.rho_sq == S(1):
                    self.rho_sq = None
            else:
                raise ValueError('In Spinor input for Lt, S*S.rev() not a scalar!\n')

        elif isinstance(mat_rep, str):  # String input
            Amat = Symbolic_Matrix(mat_rep, coords=self.Ga.coords, mode=self.mode, f=self.fct_flg)
            self.__init__(Amat, ga=self.Ga)

        else:  # Linear multivector function input
            # F is a multivector function to be tested for linearity
            F = mat_rep
            a = mv.Mv('a', 'vector', ga=self.Ga)
            b = mv.Mv('b', 'vector', ga=self.Ga)
            if F(a + b) == F(a) + F(b):
                self.lt_dict = {}
                for base in self.Ga.basis:
                    self.lt_dict[base] = (F(mv.Mv(base, ga=self.Ga))).obj
                    if not self.lt_dict[base].is_vector():
                        raise ValueError(str(mat_rep) + ' is not supported for Lt definition\n')
            else:
                raise ValueError(str(mat_rep) + ' is not supported for Lt definition\n')

    def __call__(self, v, obj=False):
        r"""
        Returns the image of the multivector :math:`A` under the linear transformation :math:`L`.

        :math:`{{L}\lp {A} \rp }` is defined by the linearity of :math:`L`, the
        vector values :math:`{{L}\lp {{{\eb}}_{i}} \rp }`, and the definition
        :math:`{{L}\lp {{{\eb}}_{i_{1}}{\wedge}\dots{\wedge}{{\eb}}_{i_{r}}} \rp } = {{L}\lp {{{\eb}}_{i_{1}}} \rp }{\wedge}\dots{\wedge}{{L}\lp {{{\eb}}_{i_{r}}} \rp }`.
        """

        if isinstance(v, mv.Mv) and self.Ga != v.Ga:
            raise ValueError('In Lt call Lt and argument refer to different vector spaces')

        if self.spinor:
            if not isinstance(v, mv.Mv):
                v = mv.Mv(v, ga=self.Ga)
            if self.rho_sq is None:
                R_v_Rrev = self.R * v * self.Rrev
            else:
                R_v_Rrev = self.rho_sq * self.R * v * self.Rrev
            if obj:
                return R_v_Rrev.obj
            else:
                return R_v_Rrev

        if isinstance(v, mv.Mv):
            if v.is_vector():
                lt_v = v.obj.xreplace(self.lt_dict)
                if obj:
                    return lt_v
                else:
                    return mv.Mv(lt_v, ga=self.Ga)
            else:
                mv_obj = v.obj
        else:
            mv_obj = mv.Mv(v, ga=self.Ga).obj

        if self.mv_dict is None:  # Build dict for linear transformation of multivector
            self.mv_dict = copy(self.lt_dict)
            for key in self.Ga.blades[2:]:
                for blade in key:
                    index = self.Ga.indexes_to_blades_dict.inverse[blade]
                    lt_blade = self(self.Ga.basis[index[0]], obj=True)
                    for i in index[1:]:
                        lt_blade = self.Ga.wedge(lt_blade, self(self.Ga.basis[i], obj=True))
                    self.mv_dict[blade] = lt_blade

        lt_v = mv_obj.xreplace(self.mv_dict)
        if obj:
            return lt_v
        else:
            return mv.Mv(lt_v, ga=self.Ga)

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

    def det(self) -> Expr:  # det(L) defined by L(I) = det(L)I
        r"""
        Returns the determinant (a scalar) of the linear transformation,
        :math:`L`, defined by :math:`{{\det}\lp {L} \rp }I = {{L}\lp {I} \rp }`.
        """

        lt_I = self(self.Ga.i, obj=True)
        det_lt_I = lt_I.subs(self.Ga.i.obj, S(1))
        return det_lt_I

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

    def adj(self) -> 'Lt':
        r"""
        Returns the adjoint (a linear transformation) of the linear
        transformation, :math:`L`, defined by :math:`a\cdot{{L}\lp {b} \rp } = b\cdot{{\bar{L}}\lp {a} \rp }`
        where :math:`a` and :math:`b` are any two vectors in the tangent space
        and :math:`\bar{L}` is the adjoint of :math:`L`.
        """

        self_adj = []
        for e_j in self.Ga.basis:
            s = S(0)
            for e_i, er_i in zip(self.Ga.basis, self.Ga.r_basis):
                s += er_i * self.Ga.hestenes_dot(e_j, self(e_i, obj=True))
            if self.Ga.is_ortho:
                self_adj.append(expand(s))
            else:
                self_adj.append(expand(s) / self.Ga.e_sq)
        return Lt(self_adj, ga=self.Ga)

    def inv(self):
        if self.spinor:
            Lt_inv = Lt(self.Rrev, ga=self.Ga)
            Lt_inv.rho_sq = S(1)/(self.rho_sq**2)
        else:
            raise ValueError('Lt inverse currently implemented only for spinor!\n')
        return Lt_inv

    def _sympystr(self, print_obj):

        if self.spinor:
            return 'R = ' + print_obj.doprint(self.R)
        else:
            pre = 'Lt('
            s = ''
            for base in self.Ga.basis:
                if base in self.lt_dict:
                    s += pre + print_obj.doprint(base) + ') = ' + print_obj.doprint(mv.Mv(self.lt_dict[base], ga=self.Ga)) + '\n'
                else:
                    s += pre + print_obj.doprint(base) + ') = 0\n'
            return s[:-1]

    def _latex(self, print_obj):

        if self.spinor:
            s = '\\left \\{ \\begin{array}{ll} '
            for base in self.Ga.basis:
                str_base = print_obj.doprint(base)
                s += 'L \\left ( ' + str_base + '\\right ) =& ' + print_obj.doprint(self.R * mv.Mv(base, ga=self.Ga) * self.Rrev) + ' \\\\ '
            s = s[:-3] + ' \\end{array} \\right \\} \n'
            return s
        else:
            s = '\\left \\{ \\begin{array}{ll} '
            for base in self.Ga.basis:
                str_base = print_obj.doprint(base)
                if base in self.lt_dict:
                    s += 'L \\left ( ' + str_base + '\\right ) =& ' + print_obj.doprint(mv.Mv(self.lt_dict[base], ga=self.Ga)) + ' \\\\ '
                else:
                    s += 'L \\left ( ' + str_base + '\\right ) =& 0 \\\\ '
            s = s[:-3] + ' \\end{array} \\right \\} \n'
            return s

    def Fmt(self, fmt=1, title=None) -> printer.GaPrintable:
        return printer._FmtResult(self, title)

    def matrix(self) -> Matrix:
        r"""
        Returns the matrix representation of the linear transformation,
        :math:`L`, defined by :math:`{{L}\lp {{{\eb}}_{i}} \rp } = L_{ij}{{\eb}}_{j}`
        where :math:`L_{ij}` is the matrix representation.
        """
        if self.mat is not None:
            return self.mat
        else:
            if self.spinor:
                self.lt_dict = {}
                for base in self.Ga.basis:
                    self.lt_dict[base] = self(base).simplify()
                self.spinor = False
                mat = self.matrix()
                self.spinor = True
                return mat
            else:
                """
                mat_rep = []
                for base in self.Ga.basis:
                    if base in self.lt_dict:
                        row = []
                        image = (self.lt_dict[base])
                        if isinstance(image, mv.Mv):
                            image = image.obj
                        coefs, bases = metric.linear_expand(image)
                        for base in self.Ga.basis:
                            try:
                                i = bases.index(base)
                                row.append(coefs[i])
                            except:
                                row.append(0)
                        mat_rep.append(row)
                    else:
                        mat_rep.append(self.Ga.n * [0])
                return Matrix(mat_rep).transpose()
                """
                self.mat = Dictionary_to_Matrix(self.lt_dict, self.Ga) * self.Ga.g
                return self.mat


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
        return print_obj.doprint(self.fvalue)

    def _latex(self, print_obj):
        if self.nargs <= 1:
            return print_obj.doprint(self.fvalue)
        expr_lst = Mlt.expand_expr(self.fvalue, self.Ga)
        latex_str = '\\begin{aligned} '
        first = True
        lcnt = print_obj._settings['galgebra_mlt_lcnt']
        cnt = 1  # Component count on line
        for term in expr_lst:
            coef_str = str(term[0])
            coef_latex = print_obj.doprint(term[0])
            term_add_flg = isinstance(term[0], Add)
            if term_add_flg:
                coef_latex = r'\left ( ' + coef_latex + r'\right ) '
            if first:
                first = False
            else:
                if coef_str[0].strip() != '-' or term_add_flg:
                    coef_latex = ' + ' + coef_latex
            for aij in term[1]:
                coef_latex += print_obj.doprint(aij) + ' '
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
            coef = S(1)
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
            self.fvalue = S(0)
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
                args, _varargs, _kwargs, _defaults = inspect.getargspec(f)
                self.nargs = len(args)
                self.f = f
                Mlt.increment_slots(self.nargs, Ga)
                self.fvalue = f(*tuple(Ga._mlt_a[0:self.nargs]))
            else:  # Tensor defined by component expression
                self.f = None
                self.nargs = len(args)
                Mlt.increment_slots(self.nargs, Ga)
                self.fvalue = f

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
