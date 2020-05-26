"""
Metric Tensor and Derivatives of Basis Vectors.
"""

import copy
import warnings
from typing import List, Optional

from sympy import (
    diff, trigsimp, Matrix, Rational,
    sqf_list, sqrt, eye, S, expand, Mul,
    Add, simplify, Expr, Function, MatrixSymbol
)

from . import printer
from ._utils import cached_property as _cached_property
from .atoms import (
    BasisVectorSymbol, DotProductSymbol, MatrixFunction, Determinant,
)

half = Rational(1, 2)


def apply_function_list(f, x):
    if isinstance(f, (tuple, list)):
        fx = x
        for fi in f:
            fx = fi(fx)
        return fx
    else:
        return f(x)


def linear_expand(expr):
    """
    linear_expand takes an expression that is the sum of a scalar
    expression and a linear combination of noncommutative terms with
    scalar coefficients and generates lists of coefficients and
    noncommutative symbols the coefficients multiply.  The list of
    noncommutatives symbols contains the scalar 1 if there is a scalar
    term in the sum and also does not contain any repeated noncommutative
    symbols.
    """
    if not isinstance(expr, Expr):
        raise TypeError('{!r} is not a SymPy Expr'.format(expr))

    expr = expand(expr)

    if expr == 0:
        coefs = [expr]
        bases = [S(1)]
        return (coefs, bases)

    if isinstance(expr, Add):
        args = expr.args
    else:
        if expr.is_commutative:
            return ([expr], [S(1)])
        else:
            args = [expr]
    coefs = []
    bases = []
    for term in args:
        if term.is_commutative:
            if S(1) in bases:
                coefs[bases.index(S(1))] += term
            else:
                bases.append(S(1))
                coefs.append(term)
        else:
            c, nc = term.args_cnc()
            base = nc[0]
            coef = Mul._from_args(c)
            if base in bases:
                coefs[bases.index(base)] += coef
            else:
                bases.append(base)
                coefs.append(coef)
    return (coefs, bases)


def linear_expand_terms(expr):
    coefs, bases = linear_expand(expr)
    return zip(coefs, bases)


def collect(A, nc_list):
    """
    Parameters
    -----------
    A :
        a linear combination of noncommutative symbols with scalar
        expressions as coefficients
    nc_list :
        noncommutative symbols in A to combine

    Returns
    -------
    sympy.Basic
        A sum of the terms containing the noncommutative symbols in `nc_list` such that no elements
        of `nc_list` appear more than once in the sum. All coefficients of a given element of `nc_list`
        are combined into a single coefficient.
    """
    coefs, bases = linear_expand(A)
    C = S(0)
    for x in nc_list:
        if x in bases:
            i = bases.index(x)
            bases.pop(i)
            C += coefs.pop(i)*x

    # add whatever is left
    for c, b in zip(coefs, bases):
        C += c * b

    return C


def square_root_of_expr(expr):
    """
    If expression is product of even powers then every power is divided
    by two and the product is returned.  If some terms in product are
    not even powers the sqrt of the absolute value of the expression is
    returned.  If the expression is a number the sqrt of the absolute
    value of the number is returned.
    """
    if expr.is_number:
        if expr > 0:
            return sqrt(expr)
        else:
            return sqrt(-expr)
    else:
        expr = trigsimp(expr)
        coef, pow_lst = sqf_list(expr)
        if coef != S(1):
            if coef.is_number:
                coef = square_root_of_expr(coef)
            else:
                coef = sqrt(abs(coef))  # Product coefficient not a number
        for p in pow_lst:
            f, n = p
            if n % 2 != 0:
                return sqrt(abs(expr))  # Product not all even powers
            else:
                coef *= f ** (n / S(2))  # Positive sqrt of the square of an expression
        return coef


def symbols_list(s, indices=None, sub=True, commutative=False):
    """
    Convert a string to a list of symbols.

    If :class:`galgebra.printer.Eprint` is enabled, the symbol names will
    contain ANSI escape sequences.

    Parameters
    ----------
    s : str
        Specification. If `indices` is specified, then this is just a prefix.
        If `indices` is not specified then this is a string of one of the forms:

        * ``prefix + "*" + index_1 + "|" + index_2 + "|" + ... + index_n``
        * ``prefix + "*" + n_indices``
        * ``name_1 + "," + name_2 + "," + ... + name_n``
        * ``name_1 + " " + name_2 + " " + ... + name_n``

    indices : list, optional
        List of indices to append to the prefix.
    sub : bool
        If true, mark as subscript separating prefix and suffix with ``_``, else
        mark as superscript using ``__``.
    commutative : bool
        Passed on to :class:`sympy.Symbol`.

    Returns
    -------
    symbols : list of :class:`sympy.Symbol`

    Examples
    --------

    Names can be comma or space separated:

    >>> symbols_list('a,b,c')
    [a, b, c]
    >>> symbols_list('a b c')
    [a, b, c]

    Mixing commas and spaces gives surprising results:

    >>> symbols_list('a b,c')
    [a b, c]

    Subscripts will be converted to superscripts if requested:

    >>> symbols_list('a_1 a_2', sub=False)
    [a__1, a__2]
    >>> symbols_list('a__1 a__2', sub=False)
    [a___1, a___2]

    But not vice versa:

    >>> symbols_list('a__1 a__2', sub=True)
    [a__1, a__2]

    Asterisk can be used for repetition:

    >>> symbols_list('a*b|c|d')
    [a_b, a_c, a_d]
    >>> symbols_list('a*3')
    [a_0, a_1, a_2]
    >>> symbols_list('a*3')
    [a_0, a_1, a_2]

    Or the indices argument:

    >>> symbols_list('a', [2, 4, 6])
    [a_2, a_4, a_6]
    >>> symbols_list('a', [2, 4, 6], sub=False)
    [a__2, a__4, a__6]

    See also
    --------
    :func:`sympy.symbols`: a similar function builtin to sympy
    """

    if isinstance(s, list):  # s is already a list of symbols
        return s

    if sub is True:  # subscripted list
        pos = '_'
    else:  # superscripted list
        pos = '__'

    if indices is None:  # symbol list completely generated by s
        if '*' in s:
            [base, index] = s.split('*')
            if '|' in s:
                index = index.split('|')
                s_lst = [base + pos + i for i in index]
            else:  # symbol list indexed with integers 0 to n-1
                try:
                    n = int(index)
                except ValueError:
                    raise ValueError(index + 'is not an integer')
                s_lst = [base + pos + str(i) for i in range(n)]
        else:
            if ',' in s:
                s_lst = s.split(',')
            else:
                s_lst = s.split(' ')
            if not sub:
                s_lst = [x.replace('_', '__', 1) for x in s_lst]

    else:  # indices symbol list used for sub/superscripts of generated symbol list
        s_lst = [s + pos + str(i) for i in indices]
    return [BasisVectorSymbol(s, commutative=commutative) for s in s_lst]


class Simp:
    modes = [simplify]

    @staticmethod
    def profile(s):
        Simp.modes = s

    @staticmethod
    def apply(expr):
        obj = S(0)
        for coef, base in linear_expand_terms(expr):
            obj += apply_function_list(Simp.modes, coef) * base
        return obj

    @staticmethod
    def applymv(mv):
        return Mv(Simp.apply(mv.obj), ga=mv.Ga)


class Metric(object):
    """
    Metric specification

    Attributes
    ----------

    g : sympy matrix[,]
        metric tensor
    g_inv : sympy matrix[,]
        inverse of metric tensor
    norm : list of sympy numbers
        normalized diagonal metric tensor
    coords : list[] of sympy symbols
        coordinate variables
    is_ortho : bool
        True if basis is orthogonal
    connect_flg : bool
        True if connection is non-zero
    basis : list[] of non-commutative sympy variables
        basis vector symbols
    r_symbols : list[] of non-commutative sympy variables
        reciprocal basis vector symbols
    n : integer
        dimension of vector space/manifold
    n_range :
        list of basis indices
    de : list[][]
        derivatives of basis functions.  Two dimensional list. First
        entry is differentiating coordiate. Second entry is basis
        vector.  Quantities are linear combinations of basis vector
        symbols.
    sig : Tuple[int, int]
        Signature of metric ``(p,q)`` where ``n = p+q``.  If metric tensor
        is numerical and orthogonal it is calculated.  Otherwise the
        following inputs are used:

        =========   ===========  ==================================
        Input       Signature     Type
        =========   ===========  ==================================
        ``"e"``     ``(n,0)``    Euclidean
        ``"m+"``    ``(n-1,1)``  Minkowski (One negative square)
        ``"m-"``    ``(1,n-1)``  Minkowski (One positive square)
        ``p``       ``(p,n-p)``  General (integer not string input)
        =========   ===========  ==================================

    gsym : str
        String for symbolic metric determinant.  If self.gsym = 'g'
        then det(g) is sympy scalar function of coordinates with
        name 'det(g)'.  Useful for complex non-orthogonal coordinate
        systems or for calculations with general metric.
    """

    count = 1

    @staticmethod
    def dot_orthogonal(V1, V2, g=None):
        """
        Returns the dot product of two vectors in an orthogonal coordinate
        system.  V1 and V2 are lists of sympy expressions.  g is
        a list of constants that gives the signature of the vector space to
        allow for non-euclidian vector spaces.

        This function is only used to form the dot product of vectors in the
        embedding space of a vector manifold or in the case where the basis
        vectors are explicitly defined by vector fields in the embedding
        space.

        A g of None is for a Euclidian embedding space.
        """
        if g is None:
            dot = 0
            for v1, v2 in zip(V1, V2):
                dot += v1 * v2
            return dot
        else:
            if len(g) == len(V1):
                dot = 0
                for v1, v2, gii in zip(V1, V2, g):
                    dot += v1 * v2 * gii
                return dot
            else:
                raise ValueError('In dot_orthogonal dimension of metric ' +
                                 'must equal dimension of vector')

    def _build_metric_element(self, s, i1, i2):
        """ Build an element for the metric of `bases[i1] . basis[i2]` """
        if s == '#':
            if i1 <= i2:  # for default element ensure symmetry
                return DotProductSymbol(self.basis[i1], self.basis[i2])
            else:
                return DotProductSymbol(self.basis[i2], self.basis[i1])
        else:  # element is fraction or integer
            return Rational(s)

    def metric_symbols_list(self, s=None):  # input metric tensor as string
        """
        rows of metric tensor are separated by "," and elements
        of each row separated by " ".  If the input is a single
        row it is assummed that the metric tensor is diagonal.

        Output is a square matrix.
        """
        if s is None:
            s = self.n * '# '
            s = self.n * (s[:-1] + ',')
            s = s[:-1]

        if isinstance(s, str):
            rows = s.split(',')
            n_rows = len(rows)

            if n_rows == 1:  # orthogonal metric
                m_lst = s.split(' ')
                m = [
                    self._build_metric_element(s, i, i)
                    for i, s in enumerate(m_lst)
                ]

                if len(m) != self.n:
                    raise ValueError('Input metric "' + s + '" has' +
                                     ' different rank than bases "' + str(self.basis) + '"')
                diagonal = eye(self.n)

                for i in self.n_range:
                    diagonal[i, i] = m[i]
                return diagonal

            else:  # non orthogonal metric
                rows = s.split(',')
                n_rows = len(rows)
                m_lst = []
                for row in rows:
                    cols = row.strip().split(' ')
                    n_cols = len(cols)
                    if n_rows != n_cols:  # non square metric
                        raise ValueError("'" + s + "' does not represent square metric")
                    m_lst.append(cols)
                n = len(m_lst)
                if n != self.n:
                    raise ValueError('Input metric "' + s + '" has' +
                                     ' different rank than bases "' + str(self.basis) + '"')
                return Matrix([
                    [
                        self._build_metric_element(s, i1, i2)
                        for i2, s in enumerate(row)
                    ]
                    for i1, row in enumerate(m_lst)
                ])

    def derivatives_of_g(self):
        # galgebra 0.5.0
        warnings.warn(
            "Metric.derivatives_of_g is deprecated, and now does nothing. "
            "the `.dg` property is now always available.",
            DeprecationWarning, stacklevel=2)

    @_cached_property
    def dg(self) -> List[List[List[Expr]]]:
        # dg[i][j][k] = \partial_{x_{k}}g_{ij}
        return [[[
            diff(self.g[i, j], x_k)
            for x_k in self.coords]
            for j in self.n_range]
            for i in self.n_range]

    @_cached_property
    def connect_flg(self) -> bool:
        """ True if connection is non-zero """
        if self.coords is None:
            return False
        else:
            return any(
                self.dg[i][j][k] != 0
                for i in self.n_range
                for j in self.n_range
                for k in self.n_range
            )

    @_cached_property
    def de(self) -> Optional[List[List[Expr]]]:
        # Derivatives of basis vectors from Christoffel symbols

        n_range = self.n_range

        if not self.connect_flg:
            return None

        # Christoffel symbols of the first kind, \Gamma_{ijk}
        # TODO handle None
        dG = self.Christoffel_symbols(mode=1)

        # de[i][j] = \partial_{x_{i}}e^{x_{j}}
        # \frac{\partial e_{j}}{\partial x^{i}} = \Gamma_{ijk} e^{k}
        de = [[
            sum([Gamma_ijk * e__k for Gamma_ijk, e__k in zip(dG[i][j], self.r_symbols)])
            for j in n_range
        ] for i in n_range]

        if self.debug:
            printer.oprint('D_{i}e^{j}', de)

        return de

    def inverse_metric(self) -> None:
        # galgebra 0.5.0
        warnings.warn(
            "Metric.inverse_metric is deprecated, and now does nothing. "
            "the `.g_inv` property is now always available.",
            DeprecationWarning, stacklevel=2)

    @_cached_property
    def g_inv(self) -> Matrix:
        """ Inverse of g """
        if self.is_ortho:  # Orthogonal metric
            g_inv = eye(self.n)
            for i in range(self.n):
                g_inv[i, i] = S(1)/self.g(i, i)
            return g_inv
        elif self.gsym is None:
            return simplify(self.g.inv())
        else:
            return self.g_adj/self.detg

    @_cached_property
    def g_adj(self) -> Matrix:
        """ Adjugate of g """
        return simplify(self.g.adjugate())

    def Christoffel_symbols(self, mode=1):
        """
        mode = 1  Christoffel symbols of the first kind
        mode = 2  Christoffel symbols of the second kind
        """

        # See if connection is zero
        if not self.connect_flg:
            return

        n_range = self.n_range

        # dg[i][j][k] = \partial_{x_{k}}g_{ij}
        dg = self.dg

        if mode == 1:

            # Christoffel symbols of the first kind, \Gamma_{ijk}
            # \partial_{x^{i}}e_{j} = \Gamma_{ijk}e^{k}

            def Gamma_ijk(i, j, k):
                return half * (dg[j][k][i] + dg[i][k][j] - dg[i][j][k])

            # dG[i][j][k] = half * (dg[j][k][i] + dg[i][k][j] - dg[i][j][k])
            dG = [[[
                Simp.apply(Gamma_ijk(i, j, k))
                for k in n_range]
                for j in n_range]
                for i in n_range]

            if self.debug:
                printer.oprint('Gamma_{ijk}', dG)
            return dG

        elif mode == 2:
            # TODO handle None
            Gamma1 = self.Christoffel_symbols(mode=1)

            # Christoffel symbols of the second kind, \Gamma_{ij}^{k} = \Gamma_{ijl}g^{lk}
            # \partial_{x^{i}}e_{j} = \Gamma_{ij}^{k}e_{k}

            def Gamma2_ijk(i, j, k):
                return sum([Gamma_ijl * self.g_inv[l, k] for l, Gamma_ijl in enumerate(Gamma1[i][j])])

            Gamma2 = [[[
                Simp.apply(Gamma2_ijk(i, j, k))
                for k in n_range]
                for j in n_range]
                for i in n_range]

            return Gamma2
        else:
            raise ValueError('In Christoffle_symobols mode = ' + str(mode) + ' is not allowed\n')

    def normalize_metric(self):

        if self.de is None:
            return

        #  Generate mapping for renormalizing reciprocal basis vectors
        renorm = [
            (self.r_symbols[ib], self.r_symbols[ib] / self.e_norm[ib])
            for ib in self.n_range  # e^{ib} --> e^{ib}/|e_{ib}|
        ]

        # Normalize derivatives of basis vectors

        for x_i in self.n_range:
            for jb in self.n_range:
                self.de[x_i][jb] = Simp.apply((((self.de[x_i][jb].subs(renorm)
                                              - diff(self.e_norm[jb], self.coords[x_i]) *
                                              self.basis[jb]) / self.e_norm[jb])))
        if self.debug:
            for x_i in self.n_range:
                for jb in self.n_range:
                    print(r'\partial_{' + str(self.coords[x_i]) + r'}\hat{e}_{' + str(self.coords[jb]) + '} =', self.de[x_i][jb])

        # Normalize metric tensor

        for ib in self.n_range:
            for jb in self.n_range:
                self.g[ib, jb] = Simp.apply(self.g[ib, jb] / (self.e_norm[ib] * self.e_norm[jb]))

        if self.debug:
            printer.oprint('e^{i}->e^{i}/|e_{i}|', renorm)
            printer.oprint('renorm(g)', self.g)

    def signature(self):
        if self.is_ortho:
            p = 0
            q = 0
            for i in self.n_range:
                g_ii = self.g[i, i]
                if g_ii.is_number:
                    if g_ii > 0:
                        p += 1
                    else:
                        q += 1
                else:
                    break
            if p + q == self.n:
                self.sig = (p, q)
                return
        if isinstance(self.sig, int):  # General signature
            if self.sig <= self.n:
                self.sig = (self.sig, self.n - self.sig)
                return
            else:
                raise ValueError('self.sig = ' + str(self.sig) + ' > self.n, not an allowed hint')
        if isinstance(self.sig, str):
            if self.sig == 'e':  # Euclidean metric signature
                self.sig = (self.n, 0)
            elif self.sig == 'm+':  # Minkowski metric signature (n-1,1)
                self.sig = (self.n - 1, 1)
            elif self.sig == 'm-':  # Minkowski metric signature (1,n-1)
                self.sig = (1, self.n - 1)
            else:
                raise ValueError('self.sig = ' + str(self.sig) + ' is not an allowed hint')
            return
        raise ValueError(str(self.sig) + ' is not allowed value for self.sig')

    @_cached_property
    def detg(self) -> Expr:
        r""" Determinant of :math:`g`, :math:`\det g` """
        if self.gsym is None:
            g = self.g
        else:
            # Define name of metric tensor determinant as sympy symbol
            if self.coords is None:
                g = MatrixSymbol(self.gsym, self.n, self.n)
            else:
                g = MatrixFunction(self.gsym, self.n, self.n)(*self.coords)
        return Determinant(g)

    def __init__(
        self, basis, *,
        g=None,
        coords=None,
        X=None,
        norm=False,
        debug=False,
        gsym=None,
        sig='e',
        Isq='-'
    ):
        """
        Parameters
        ----------
        basis :
            string specification
        g :
            metric tensor
        coords :
            manifold/vector space coordinate list/tuple  (sympy symbols)
        X :
            vector manifold function
        norm :
            True to normalize basis vectors
        debug :
            True to print out debugging information
        gsym :
            String s to use ``"det("+s+")"`` function in reciprocal basis
        sig :
            Signature of metric, default is (n,0) a Euclidean metric
        Isq :
            Sign of square of pseudo-scalar, default is '-'
        """

        self.name = 'GA' + str(Metric.count)
        Metric.count += 1

        if not isinstance(basis, str):
            raise TypeError('"' + str(basis) + '" must be string')

        self.sig = sig  # Hint for metric signature
        self.gsym = gsym
        self.Isq = Isq  #: Sign of I**2, only needed if I**2 not a number

        self.debug = debug
        self.is_ortho = False  # Is basis othogonal
        self.coords = coords  # Manifold coordinates
        self.norm = norm  # True to normalize basis vectors
        # Generate list of basis vectors and reciprocal basis vectors
        # as non-commutative symbols

        if ' ' in basis or ',' in basis or '*' in basis:  # bases defined by substrings separated by spaces or commas
            self.basis = symbols_list(basis)
            self.r_symbols = symbols_list(basis, sub=False)
        else:
            if coords is not None:  # basis defined by root string with symbol list as indices
                self.basis = symbols_list(basis, coords)
                self.r_symbols = symbols_list(basis, coords, sub=False)
                self.coords = coords
                if self.debug:
                    printer.oprint('x^{i}', self.coords)
            else:
                raise ValueError('for basis "' + basis + '" coords must be entered')

        if self.debug:
            printer.oprint('e_{i}', self.basis, 'e^{i}', self.r_symbols)
        self.n = len(self.basis)
        self.n_range = list(range(self.n))

        # Generate metric as list of lists of symbols, rationals, or functions of coordinates

        if g is None:
            if X is None:  # default metric from dot product of basis as symbols
                self.g = self.metric_symbols_list()
            else:  # Vector manifold
                if coords is None:
                    raise ValueError('For metric derived from vector field ' +
                                     ' coordinates must be defined.')
                else:  # Vector manifold defined by vector field
                    # Get basis vectors by differentiating vector field
                    dX = [
                        [diff(x, coord) for x in X]
                        for coord in coords
                    ]
                    self.g = Matrix([
                        [
                            trigsimp(Metric.dot_orthogonal(dx1, dx2, g))
                            for dx2 in dX
                        ]
                        for dx1 in dX
                    ])
                    if self.debug:
                        printer.oprint('X_{i}', X, 'D_{i}X_{j}', dX)

        else:  # metric is symbolic or list of lists of functions of coordinates
            if isinstance(g, str):  # metric elements are symbols or constants
                if g == 'g':  # general symbolic metric tensor (g_ij functions of position)
                    self.g = Matrix([
                        [
                            Function('g_{}_{}'.format(coord, coord2))(*self.coords)
                            for coord2 in self.coords
                        ]
                        for coord in self.coords
                    ])
                    self.g_inv = Matrix([
                        [
                            Function('g__{}__{}'.format(coord, coord2))(*self.coords)
                            for coord2 in self.coords
                        ]
                        for coord in self.coords
                    ])
                else:  # specific symbolic metric tensor (g_ij are symbolic or numerical constants)
                    self.g = self.metric_symbols_list(g)  # construct symbolic metric from string and basis
            else:  # metric is given as list of function or list of lists of function or matrix of functions
                if isinstance(g, Matrix):
                    self.g = g
                else:
                    if isinstance(g[0], list):
                        self.g = Matrix(g)
                    else:
                        m = eye(len(g))
                        for i in range(len(g)):
                            m[i, i] = g[i]
                        self.g = m

        self.g_raw = copy.deepcopy(self.g)  # save original metric tensor for use with submanifolds

        if self.debug:
            printer.oprint('g', self.g)

        # Determine if metric is orthogonal

        self.is_ortho = all(
            self.g[i, j] == 0
            for i in self.n_range
            for j in self.n_range
            if i < j
        )

        self.g_is_numeric = all(
            self.g[i, j].is_number
            for i in self.n_range
            for j in self.n_range
            if i < j
        )

        if self.coords is not None:
            if self.norm:  # normalize basis, metric, and derivatives of normalized basis
                if not self.is_ortho:
                    raise ValueError('!!!!Basis normalization only implemented for orthogonal basis!!!!')
                self.e_norm = [
                    square_root_of_expr(self.g[i, i])
                    for i in self.n_range
                ]
                if debug:
                    printer.oprint('|e_{i}|', self.e_norm)
            else:
                self.e_norm = None

        if self.norm:
            if self.is_ortho:
                self.normalize_metric()
            else:
                raise ValueError('!!!!Basis normalization only implemented for orthogonal basis!!!!')

        if not self.g_is_numeric:
            self.signature()
            # Sign of square of pseudo scalar
            self.e_sq_sgn = '+'
            if ((self.n*(self.n-1))//2+self.sig[1]) % 2 == 1:
                self.e_sq_sgn = '-'

        if self.debug:
            print('signature =', self.sig)
