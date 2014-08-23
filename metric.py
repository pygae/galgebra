#metric.py

import sys
import copy
import itertools
from sympy import diff, trigsimp, Matrix, Rational, \
    sqf_list, Symbol, sqrt, eye, S, expand, Mul, \
    Add, simplify, together, ratsimp, Expr, latex, \
    numbers

import printer

half = Rational(1, 2)

def in_ipynb():
    try:
        cfg = get_ipython().config
        if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
            return True
        else:
            return False
    except NameError:
        return False

def str_to_lst(s):
    if '[' in s:
        s = s.replace('[', '')
    if ']' in s:
        s = s.replace(']', '')
    s_lst = s.split(',')
    v_lst = []
    for x in s_lst:
        try:
            v_lst.append(int(s))
        except ValueError:
            v_lst.append(Symbol(s, real=True))
    return v_lst


def linear_expand(expr, mode=True):

    if isinstance(expr, Expr):
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
    if mode:
        return (coefs, bases)
    else:
        return zip(coefs, bases)


def square_root_of_expr(expr):
    if expr.is_number:
        if expr > 0:
            return(sqrt(expr))
        else:
            return(sqrt(-expr))
    else:
        expr = trigsimp(expr)
        (coef, pow_lst) = sqf_list(expr)
        if coef != S(1):
            if coef.is_number:
                coef = square_root_of_expr(coef)
            else:
                coef = sqrt(abs(coef))
        for p in pow_lst:
            (f, n) = p
            if n % 2 != 0:
                return(sqrt(abs(expr)))
            else:
                coef *= f ** (n / 2)
        return coef


def symbols_list(s, indices=None, sub=True, commutative=False):

    if isinstance(s, list):  # s is already a list of symbols
        return(s)

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
    return [Symbol(printer.Eprint.Base(s), commutative=commutative) for s in s_lst]


def test_init_slots(init_slots, **kwargs):
    """
    Tests kwargs for allowed keyword arguments as defined by dictionary
    init_slots.  If keyword argument defined by init_slots is not present
    set default value asdefined by init_slots.  Allow for backward
    compatible keyword arguments by equivalencing keywords by setting
    default value of backward compatible keyword to new keyword and then
    referencing new keywork (see init_slots for Metric class and equivalence
    between keywords 'g' and 'metric')
    """

    for slot in kwargs:
        if slot not in init_slots:
            print 'Allowed keyed input arguments'
            for key in init_slots:
                print key + ': ' + init_slots[key][1]
            raise ValueError('"' + slot + ' = " not in allowed values.')
    for slot in init_slots:
        if slot in kwargs:
            if init_slots[slot][0] in init_slots:  # redirect for backward compatibility
                kwargs[init_slots[slot][0]] = kwargs[slot]
        else:  # use default value
            if init_slots[slot][0] in init_slots:  # redirect for backward compatibility
                kwargs[init_slots[slot][0]] = init_slots[init_slots[slot][0]][0]
            kwargs[slot] = init_slots[slot][0]
    return kwargs


class Simp:
    modes = (simplify, trigsimp)

    @staticmethod
    def profile(s):
        Simp.modes = s
        return

    @staticmethod
    def apply(expr):
        (coefs, bases) = linear_expand(expr)
        obj = S(0)
        if isinstance(Simp.modes, list) or isinstance(Simp.modes, tuple):
            for (coef, base) in zip(coefs, bases):
                for mode in Simp.modes:
                    coef = mode(coef)
                obj += coef * base
        else:
            for (coef, base) in zip(coefs, bases):
                obj += modes(coef) * base
        return obj

    @staticmethod
    def applymv(mv):
        (coefs, bases) = linear_expand(mv.obj)
        obj = S(0)
        if isinstance(Simp.modes, list) or isinstance(Simp.modes, tuple):
            for (coef, base) in zip(coefs, bases):
                for mode in Simp.modes:
                    coef = mode(coef)
                obj += coef * base
        else:
            for (coef, base) in zip(coefs, bases):
                obj += modes(coef) * base
        mv.obj = obj
        return mv


class Metric(object):

    count = 1

    init_slots = {'g': (None, 'metric tensor'),
                  'coords': (None, 'manifold/vector space coordinate list/tuple'),
                  'X': (None, 'vector manifold function'),
                  'norm': (False, 'True to normalize basis vectors'),
                  'debug': (False, 'True to print out debugging information'),
                  'gsym': (None, 'String s to use "det("+s+")" function in reciprocal basis')}

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
            for (v1, v2) in zip(V1, V2):
                dot += v1 * v2
            return dot
        else:
            if len(g) == len(V1):
                dot = 0
                for (v1, v2, gii) in zip(V1, V2, g):
                    dot += v1 * v2 * gii
                return dot
            else:
                raise ValueError('In dot_orthogonal dimension of metric ' +
                                 'must equal dimension of vector')

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

        if isinstance(s, basestring):
            rows = s.split(',')
            n_rows = len(rows)

            if n_rows == 1:  # orthogonal metric
                m_lst = s.split(' ')
                m = []
                for (s, base) in zip(m_lst, self.basis):
                    if s == '#':
                        s_symbol = Symbol('(' + str(base) + '.' + str(base) + ')', real=True)
                    else:
                        if '/' in s:
                            [num, dem] = s.split('/')
                            s_symbol = Rational(num, dem)
                        else:
                            s_symbol = Rational(s)
                    m.append(s_symbol)

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
                m = []
                n = len(m_lst)
                if n != self.n:
                    raise ValueError('Input metric "' + s + '" has' +
                                     ' different rank than bases "' + str(self.basis) + '"')
                n_range = range(n)
                for (row, i1) in zip(m_lst, n_range):
                    row_symbols = []
                    for (s, i2) in zip(row, n_range):
                        if s == '#':
                            if i1 <= i2:  # for default elment insure symmetry
                                row_symbols.append(Symbol('(' + str(self.basis[i1]) +
                                                          '.' + str(self.basis[i2]) + ')', real=True))
                            else:
                                row_symbols.append(Symbol('(' + str(self.basis[i2]) +
                                                          '.' + str(self.basis[i1]) + ')', real=True))
                        else:
                            if '/' in s:  # element is fraction
                                [num, dem] = s.split('/')
                                row_symbols.append(Rational(num, dem))
                            else:  # element is integer
                                row_symbols.append(Rational(s))
                    m.append(row_symbols)
                m = Matrix(m)
                return m

    def derivatives_of_basis(self):  # Derivatives of basis vectors from Christ-Awful symbols

        dg = []  # dg[i][j][k] = \partial_{x_{k}}g_{ij}

        for i in self.n_range:
            dg_row = []
            for j in self.n_range:
                dg_row.append([diff(self.g[i, j], coord) for coord in self.coords])
            dg.append(dg_row)

        # See if metric is flat

        self.connect_flg = False

        for i in self.n_range:
            for j in self.n_range:
                for k in self.n_range:
                    if dg[i][j][k] != 0:
                        self.connect_flg = True
                        break

        if not self.connect_flg:
            self.de = None
            return

        n_range = range(len(self.basis))

        de = []  # de[i][j] = \partial_{x_{i}}e^{x_{j}}

        # Christoffel symbols of the first kind, \Gamma_{ijk}

        for i in n_range:
            de_row = []
            for j in n_range:
                Gamma = []
                for k in n_range:
                    gamma = half * (dg[j][k][i] + dg[i][k][j] - dg[i][j][k])
                    Gamma.append(Simp.apply(gamma))
                de_row.append(sum([gamma * base for (gamma, base) in zip(Gamma, self.r_symbols)]))
            de.append(de_row)
        if self.debug:
            printer.oprint('D_{i}e^{j}', de)
        self.de = de
        return

    def normalize_metric(self):

        if self.de is None:
            return

        renorm = []
        #  Generate mapping for renormalizing reciprocal basis vectors
        for ib in self.n_range:  # e^{ib} --> e^{ib}/|e_{ib}|
            renorm.append((self.r_symbols[ib], self.r_symbols[ib] / self.e_norm[ib]))

        # Normalize derivatives of basis vectors

        for x_i in self.n_range:
            for jb in self.n_range:
                self.de[x_i][jb] = Simp.apply((((self.de[x_i][jb].subs(renorm)
                                              - diff(self.e_norm[jb], self.coords[x_i]) *
                                              self.basis[jb]) / self.e_norm[jb])))
        if self.debug:
            for x_i in self.n_range:
                for jb in self.n_range:
                    print '\partial_{' + str(self.coords[x_i]) + '}\hat{e}_{' + str(self.coords[jb]) + '} =', self.de[x_i][jb]

        # Normalize metric tensor

        for ib in self.n_range:
            for jb in self.n_range:
                self.g[ib, jb] = Simp.apply(self.g[ib, jb] / (self.e_norm[ib] * self.e_norm[jb]))

        if self.debug:
            printer.oprint('e^{i}->e^{i}/|e_{i}|', renorm)
            printer.oprint('renorm(g)', self.g)

        return

    def __init__(self, basis, **kwargs):

        kwargs = test_init_slots(Metric.init_slots, **kwargs)

        self.name = 'GA' + str(Metric.count)
        Metric.count += 1

        if not isinstance(basis, basestring):
            raise TypeError('"' + str(basis) + '" must be string')

        X = kwargs['X']
        g = kwargs['g']
        debug = kwargs['debug']
        coords = kwargs['coords']
        norm = kwargs['norm']
        self.gsym = kwargs['gsym']

        """
        Normalization for reciprocal vectors if you do not wish to
        explicitly calculate the determinate of the metric tensor.
        """

        self.debug = debug
        self.is_ortho = False  # Is basis othogonal
        self.coords = coords  # Manifold coordinates
        if self.coords is None:
            self.connect_flg = False
        else:
            self.connect_flg = True  # Connection needed for postion dependent metric
        self.norm = norm  # True to normalize basis vectors

        # Generate list of basis vectors and reciprocal basis vectors
        # as non-commutative symbols

        if ' ' in basis or '*' in basis:  # bases defined by substrings separated by spaces
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
        self.n_range = range(self.n)

        # Generate metric as list of lists of symbols, rationals, or functions of coordinates

        if g is None:
            if X is None:  # default metric from dot product of basis as symbols
                self.g = self.metric_symbols_list()
            else:  # Vector manifold
                if coords is None:
                    raise ValueError('For metric derived from vector field ' +
                                     ' coordinates must be defined.')
                else:  # Vector manifold defined by vector field
                    dX = []
                    for coord in coords:  # Get basis vectors by differentiating vector field
                        dX.append([diff(x, coord) for x in X])
                    g_tmp = []
                    for dx1 in dX:
                        g_row = []
                        for dx2 in dX:
                            dx1_dot_dx2 = trigsimp(Metric.dot_orthogonal(dx1, dx2, g))
                            g_row.append(dx1_dot_dx2)
                        g_tmp.append(g_row)
                    self.g = Matrix(g_tmp)
                    if self.debug:
                        printer.oprint('X_{i}', X, 'D_{i}X_{j}', dX)

        else:  # metric is symbolic or list of lists of functions of coordinates
            if isinstance(g, basestring):  # metric elements are symbols or constants
                if g == 'g':  # general symbolic metric tensor (g_ij functions of position)
                    g_lst = []
                    g_inv_lst = []
                    for coord1 in self.coords:
                        i1 = str(coord2)
                        tmp = []
                        tmp_inv = []
                        for coord2 in self.coords:
                            i2 = str(coord2)
                            tmp.append(Function('g_'+i1+'_'+i2)(*self.coords))
                            tmp_inv.append(Function('g__'+i1+'__'+i2)(*self.coords))
                        g_lst.append(tmp)
                        g_inv_lst.append(tmp_inv)
                    self.g = Matrix(g_lst)
                    self.g_inv = Matrix(g_inv_lst)
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

        self.g_raw = copy.copy(self.g)  # save original metric tensor for use with submanifolds

        if self.debug:
            printer.oprint('g', self.g)

        if self.coords is not None:
            self.derivatives_of_basis()  # calculate derivatives of basis
            if self.norm:  # normalize basis, metric, and derivatives of normalized basis
                self.e_norm = []
                for i in self.n_range:
                    self.e_norm.append(square_root_of_expr(self.g[i, i]))
                if debug:
                    printer.oprint('|e_{i}|', self.e_norm)
            else:
                self.e_norm = None

        self.is_ortho = True

        # Determine if metric is orthogonal

        for i in self.n_range:
            for j in self.n_range:
                if i < j:
                    if self.g[i, j] != 0:
                        self.is_ortho = False
                        break

        if self.norm:
            self.normalize_metric()

if __name__ == "__main__":
    pass

