#mv.py

import itertools
import copy
import numbers
import operator
import functools
#from compiler.ast import flatten
from operator import itemgetter, mul, add
from itertools import combinations
from sympy import Symbol, Function, S, expand, Add, Mul, Pow, Basic, \
    sin, cos, sinh, cosh, sqrt, trigsimp, expand, \
    simplify, diff, Rational, Expr, Abs, collect, combsimp
from sympy import exp as sympy_exp
from sympy import N as Nsympy
import printer
import metric
import ga
import sys
import mv
import lt
from functools import reduce


ONE = S(1)
ZERO = S(0)

half = Rational(1, 2)

modules = \
"""
from sympy import symbols, sin, Function
from mv import Mv, pDop, sDop, Dop
from ga import Ga, half
from printer import Eprint, xdvi
from lt import Lt
"""

########################### Multivector Class ##########################


class Mv(object):
    """
    Wrapper class for multivector objects (self.obj) so that it is easy
    to overload operators (*,^,|,<,>) for the various multivector
    products and for printing.  Also provides an __init__ fuction to
    easily instanciate multivector objects.  Additionally, the functionality
    of the multivector derivative have been added via the special vector
    'grad' so that one can take the geometric derivative of a multivector
    function 'A' by applying 'grad' from the left, 'grad*A', or the
    right 'A*grad' for both the left and right derivatives.  The operator
    between the 'grad' and the 'A' can be any of the multivector product
    operators.

    If 'f' is a scalar function 'grad*f' is the usual gradient of a function.
    If 'A' is a vector function 'grad|f' is the divergence of 'A' and
    '-I*(grad^A)' is the curl of 'A' (I is the pseudo scalar for the geometric
    algebra)

    Data Variables -



    """

    ################### Multivector initialization #####################

    fmt = 1
    latex_flg = False
    restore = False
    init_slots = {'f': (False, 'True if function of coordinates'),
                  'ga': (None, 'Geometric algebra to be used with multivectors'),
                  'coords': (None, 'Coordinates to be used with multivector function'),
                  'recp': (None, 'Normalization for reciprocal vector')}
    dual_mode_lst = ['+I','I+','+Iinv','Iinv+','-I','I-','-Iinv','Iinv-']

    @staticmethod
    def setup(ga):
        """
        Set up constant mutilvectors reqired for multivector class for
        a given geometric algebra, 'ga'.
        """
        Mv.fmt = 1

        basis = [Mv(x, ga=ga) for x in ga.basis]
        I = Mv(ga.iobj, ga=ga)  # default pseudoscalar
        x = Mv('XxXx', 'vector', ga=ga)  # testing vectors
        # return default basis vectors and grad vector if coords defined
        return I, basis, x

    @staticmethod
    def get_Ga(name):
        return(Mv.Ga[name])

    @staticmethod
    def Format(mode=1):
        Mv.latex_flg = True
        Mv.fmt = mode
        return

    @staticmethod
    def Mul(A, B, op):
        """
        Function for all types of geometric multiplications called by
        overloaded operators for *, ^, |, <, and >.
        """
        if not isinstance(A, Mv):
            A = B.Ga.mv(A)
        if not isinstance(B, Mv):
            B = A.Ga.mv(B)

        if op == '*':
            return A * B
        elif op == '^':
            return A ^ B
        elif op == '|':
            return A | B
        elif op == '<':
            return A < B
        elif op == '>':
            return A > B
        else:
            raise ValeError('Operation ' + op + 'not allowed in Mv.Mul!')
        return

    def characterise_Mv(self):
        if self.char_Mv:
            return
        obj = expand(self.obj)
        if isinstance(obj, numbers.Number):
            self.i_grade = 0
            self.is_blade_rep = True
            self.grades = [0]
            return
        if  obj.is_commutative:
            self.i_grade = 0
            self.is_blade_rep = True
            self.grades = [0]
            return
        if isinstance(obj, Add):
            args = obj.args
        else:
            if obj in self.Ga.blades_lst:
                self.is_blade_rep = True
                self.i_grade = self.Ga.blades_to_grades_dict[obj]
                self.grades = [self.i_grade]
                self.char_Mv = True
                self.blade_flg = True
                return
            else:
                args = [obj]

        grades = []
        #print 'args =', args
        self.is_blade_rep = True
        for term in args:
            if term.is_commutative:
                if 0 not in grades:
                    grades.append(0)
            else:
                c, nc = term.args_cnc(split_1=False)
                blade = nc[0]
                #print 'blade =',blade
                if blade in self.Ga.blades_lst:
                    grade = self.Ga.blades_to_grades_dict[blade]
                    if not grade in grades:
                        grades.append(grade)
                else:
                    self.char_Mv = True
                    self.is_blade_rep = False
                    self.i_grade = None
                    return
        if len(grades) == 1:
            self.i_grade = grades[0]
        else:
            self.i_grade = None
        self.grades = grades
        self.char_Mv = True
        return

    def make_grade(self, *kargs, **kwargs):
    # Called by __init__ to make a pure grade multivector.

        grade = kargs[1]
        self.i_grade = grade
        if isinstance(kargs[0],str):
            root = kargs[0] + '__'
            if isinstance(kwargs['f'], bool) and not kwargs['f']:  #Is a constant mulitvector function
                self.obj = sum([Symbol(root + super_script, real=True) * base
                                for (super_script, base) in zip(self.Ga.blade_super_scripts[grade], self.Ga.blades[grade])])

            else:
                if isinstance(kwargs['f'], bool):  #Is a multivector function of all coordinates
                    self.obj = sum([Function(root + super_script, real=True)(*self.Ga.coords) * base
                        for (super_script, base) in zip(self.Ga.blade_super_scripts[grade], self.Ga.blades[grade])])
                else: #Is a multivector function of tuple kwargs['f'] variables
                    self.obj = sum([Function(root + super_script, real=True)(*kwargs['f']) * base
                        for (super_script, base) in zip(self.Ga.blade_super_scripts[grade], self.Ga.blades[grade])])
        else:
            if isinstance(kargs[0],(list,tuple)):
                if len(kargs[0]) <= len(self.Ga.blades[grade]):
                    self.obj = sum([coef * base
                        for (coef, base) in zip(kargs[0], self.Ga.blades[grade][:len(kargs[0])])])
                else:
                    pass
            else:
                pass
        return

    def make_scalar(self, *kargs, **kwargs):
    # Called by __init__ to make a scalar multivector

        if isinstance(kargs[0],str):
            if 'f' in kwargs and isinstance(kwargs['f'],bool):
                if kwargs['f']:
                    self.obj = Function(kargs[0])(*self.Ga.coords)
                else:
                    self.obj = Symbol(kargs[0], real=True)
            else:
                if 'f' in kwargs and isinstance(kwargs['f'],tuple):
                    self.obj = Function(kargs[0])(*kwargs['f'])
        else:
            self.obj = kargs[0]
        return

    def make_vector(self, *kargs, **kwargs):
    # Called by __init__ to make a vector multivector

        self.make_grade(*(kargs[0], 1), **kwargs)
        return

    def make_bivector(self, *kargs, **kwargs):
    # Called by __init__ to make a bivector multivector

        self.make_grade(*(kargs[0], 2), **kwargs)
        return

    def make_pseudo_scalar(self, *kargs, **kwargs):
    # Called by __init__ to make a pseudo scalar multivector

        self.make_grade(*(kargs[0], self.Ga.n), **kwargs)
        return

    def make_multivector(self, *kargs, **kwargs):
    # Called by __init__ to make a general (2**n components) multivector

        self.make_scalar(kargs[0], **kwargs)
        tmp = self.obj
        for grade in self.Ga.n_range:
            self.make_grade(*(kargs[0], grade + 1), **kwargs)
            tmp += self.obj
        self.obj = tmp
        return

    def make_spinor(self, *kargs, **kwargs):
    # Called by __init__ to make a general even (spinor) multivector

        self.make_scalar(kargs[0], **kwargs)
        tmp = self.obj
        for grade in self.Ga.n_range:
            if (grade + 1) % 2 == 0:
                self.make_grade(*(kargs[0], grade + 1), **kwargs)
                tmp += self.obj
        self.obj = tmp
        return

    def make_odd(self, *kargs, **kwargs):
    # Called by __init__ to make a general odd multivector
        self.make_scalar(kargs[0], **kwargs)
        tmp = S(0)
        for grade in self.Ga.n_range:
            if (grade + 1) % 2 == 1:
                self.make_grade(*(kargs[0], grade + 1), **kwargs)
                tmp += self.obj
        self.obj = tmp
        return

    init_dict = {'scalar': make_scalar,
                 'vector': make_vector,
                 'bivector': make_bivector,
                 'grade2': make_bivector,
                 'pseudo': make_pseudo_scalar,
                 'mv': make_multivector,
                 'spinor': make_spinor,
                 'even': make_spinor,
                 'odd': make_odd,
                 'grade': make_grade}

    def __init__(self, *kargs, **kwargs):

        if 'ga' not in kwargs:
            raise ValueError("Geometric algebra key inplut 'ga' required")

        kwargs = metric.test_init_slots(Mv.init_slots, **kwargs)

        self.Ga = kwargs['ga']
        self.recp = kwargs['recp']  # Normalization for reciprocal vectors

        self.char_Mv = False
        self.i_grade = None  # if pure grade mv, grade value
        self.grades = None  # list of grades in mv
        self.is_blade_rep = True  # flag for blade representation
        self.blade_flg = None  # if is_blade is called flag is set
        self.versor_flg = None  # if is_versor is called flag is set
        self.coords = self.Ga.coords
        self.title = None
        self.dot_flg = False  # Flag for over dot notation

        if len(kargs) == 0:  # default constructor 0
            self.obj = S(0)
            self.i_grade = 0
        elif len(kargs) == 1 and not isinstance(kargs[0], str):  # copy constructor
            x = kargs[0]
            if isinstance(x, Mv):
                self.obj = x.obj
                self.is_blade_rep = x.is_blade_rep
                self.i_grade = x.i_grade
            else:
                if isinstance(x, Expr):  #copy constructor for obj expression
                    self.obj = x
                else:  #copy constructor for scalar obj expression
                    self.obj = S(x)
                self.is_blade_rep = True
                self.characterise_Mv()
        else:
            if not isinstance(kargs[1],int):
                if isinstance(kargs[1],str) and kargs[1] not in Mv.init_dict:
                    raise ValueError('"' + str(kargs[1]) + '" not an allowed multivector type.')

            if isinstance(kargs[1],str):
                mode = kargs[1]
                kargs = [kargs[0]] + list(kargs[2:])
                Mv.init_dict[mode](self, *kargs, **kwargs)
            else:  # kargs[1] = r (integer) Construct grade r multivector
                if kargs[1] == 0:
                    Mv.init_dict['scalar'](self, *kargs, **kwargs)
                else:
                    Mv.init_dict['grade'](self, *kargs, **kwargs)

            #if isinstance(kargs[0],str):
            #    self.title = kargs[0]
            self.characterise_Mv()

    ################# Multivector member functions #####################

    def reflect_in_blade(self, blade):  # Reflect mv in blade
        # See Mv class functions documentation
        if blade.is_blade():
            self.characterise_Mv()
            blade.characterise_Mv()
            blade_inv = blade.rev() / blade.norm2()
            grade_dict = self.Ga.grade_decomposition(self)
            blade_grade = blade.i_grade
            reflect = Mv(0,'scalar',ga=self.Ga)
            for grade in list(grade_dict.keys()):
                if (grade * (blade_grade + 1)) % 2 == 0:
                    reflect += blade * grade_dict[grade] * blade_inv
                else:
                    reflect -= blade * grade_dict[grade] * blade_inv
            return reflect
        else:
            raise ValueError(str(blade) + 'is not a blade in reflect_in_blade(self, blade)')

    def project_in_blade(self,blade):
        # See Mv class functions documentation
        if blade.is_blade():
            blade.characterise_Mv()
            blade_inv = blade.rev() / blade.norm2()
            return (self < blade) * blade_inv  # < is left contraction
        else:
            raise ValueError(str(blade) + 'is not a blade in project_in_blade(self, blade)')

    def rotate_multivector(self,itheta,hint='-'):
        Rm = (-itheta/S(2)).exp(hint)
        Rp = (itheta/S(2)).exp(hint)
        return Rm * self * Rp

    def base_rep(self):
        if self.is_blade_rep:
            self.obj = self.Ga.blade_to_base_rep(self.obj)
            self.is_blade_rep = False
            return self
        else:
            return self

    def blade_rep(self):
        if self.is_blade_rep:
            return self
        else:
            self.obj = self.Ga.base_to_blade_rep(self.obj)
            self.is_blade_rep = True
            return self

    def __ne__(self, A):
        if isinstance(A, Mv):
            diff = (self - A).expand()
            if diff.obj == S(0):
                return False
            else:
                return True
        else:
            if self.is_scalar() and self.obj == A:
                return False
            else:
                return True

    def __eq__(self, A):
        if isinstance(A, Mv):
            diff = self - A
            if not isinstance(diff, Mv):
                if diff == 0:
                    return True
                else:
                    return False
            else:
                diff = (self - A).expand().simplify()
            if diff.obj == S(0):
                return True
            else:
                return False
        else:
            if self.is_scalar() and self.obj == A:
                return True
            else:
                return False


    """
    def __eq__(self, A):
        if not isinstance(A, Mv):
            if not self.is_scalar():
                return False
            if expand(self.obj) == expand(A):
                return True
            else:
                return False
        if self.is_blade_rep != A.is_blade_rep:
            self = self.blade_rep()
            A = A.blade_rep()
        coefs, bases = metric.linear_expand(self.obj)
        Acoefs, Abases = metric.linear_expand(A.obj)
        if len(bases) != len(Abases):
            return False
        if set(bases) != set(Abases):
            return False
        for base in bases:
            index = bases.index(base)
            indexA = Abases.index(base)
            if expand(coefs[index]) != expand(Acoefs[index]):
                return False
        return True
    """

    def __neg__(self):
        return Mv(-self.obj, ga=self.Ga)

    def __add__(self, A):

        if not isinstance(A, (Mv, Dop, Sdop, Pdop)):
            mv = Mv(self.obj + A, ga=self.Ga)
            if mv.is_zero():
                return 0
            return mv

        if self.Ga.name != A.Ga.name:
            raise ValueError('In + operation Mv arguments are not from same geometric algebra')

        if isinstance(A, Dop):
            return Dop.Add(A, self)

        if isinstance(A, Sdop):
            if self.is_scalar():
                return Sdop.Add(self.obj,A)
            return Dop.Add(A.promote(), self)

        if isinstance(A, Pdop):
            if self.is_scalar():
                return Sdop.Add(self.obj,A.promote())
            return Dop.Add(A.promote().promote(), self)

        if self.is_blade_rep == A.is_blade_rep:
            mv = Mv(self.obj + A.obj, ga=self.Ga)
            if mv.is_zero():
                return 0
            return mv
        else:
            if self.is_blade_rep:
                A = A.blade_rep()
            else:
                self = self.blade_rep()
            mv = Mv(self.obj + A.obj, ga=self.Ga)
            if mv.is_zero():
                return 0
            return mv

    def __radd__(self, A):
        return(self + A)

    def __add_ab__(self, A):  # self += A
        self.obj += A.obj
        self.char_Mv = False
        self.characterise_Mv()
        return(self)

    def __sub__(self, A):

        if isinstance(A, Pdop):
            A = A.promote()

        if isinstance(A, Sdop):
            A = A.promote()

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            mv = Mv(self.obj - A, ga=self.Ga)
            if mv.is_zero():
                return 0
            return mv

        if self.Ga.name != A.Ga.name:
            raise ValueError('In - operation Mv arguments are not from same geometric algebra')

        if isinstance(A, Dop):
            return Dop.Add(self, -A)

        if self.is_blade_rep == A.is_blade_rep:
            mv = Mv(self.obj - A.obj, ga=self.Ga)
            if mv.is_zero():
                return 0
            return mv
        else:
            if self.is_blade_rep:
                A = A.blade_rep()
            else:
                self = self.blade_rep()
            mv = Mv(self.obj - A.obj, ga=self.Ga)
            if mv.is_zero():
                return 0
            return mv

    def __rsub__(self, A):
        return -self + A

    def __sub_ab__(self, A):  # self -= A
        self.obj -= A.obj
        self.char_Mv = False
        self.characterise_Mv()
        return(self)

    def __mul__(self, A):  # Mv geometric (*) product

        if not isinstance(A, (Mv, Dop, Sdop, Pdop)):
            return Mv(expand(A * self.obj), ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In * operation Mv arguments are not from same geometric algebra')

        if isinstance(A, (Dop, Sdop, Pdop)):
            if isinstance(A, Pdop):
                A = A.promote().promote()
            if isinstance(A, Sdop):
                A = A.promote()
            #return A.Mul(self, A, op='*')
            return Dop.__rmul__(A, self)

        if self.is_scalar():
            return Mv(self.obj * A, ga=self.Ga)

        if self.is_blade_rep and A.is_blade_rep:
            self = self.base_rep()
            A = A.base_rep()

            selfxA = Mv(self.Ga.mul(self.obj, A.obj), ga=self.Ga)
            selfxA.is_blade_rep = False
            selfxA = selfxA.blade_rep()

            self = self.blade_rep()
            A = A.blade_rep()
        elif self.is_blade_rep:
            self = self.base_rep()

            selfxA = Mv(self.Ga.mul(self.obj, A.obj), ga=self.Ga)
            selfxA.is_blade_rep = False
            selfxA = selfxA.blade_rep()

            self = self.blade_rep()
        elif A.is_blade_rep:
            A = A.base_rep()

            selfxA = Mv(self.Ga.mul(self.obj, A.obj), ga=self.Ga)
            selfxA.is_blade_rep = False
            selfxA = selfxA.blade_rep()

            A = A.blade_rep()
        else:
            selfxA = Mv(self.Ga.mul(self.obj, A.obj), ga=self.Ga)

        return selfxA

    def __rmul__(self, A):  # Mv geometric (*) product
        return Mv(expand(A * self.obj), ga=self.Ga)

    def __and__(self, A):
        return((self*A).scalar())

    def __mul_ab__(self, A):  # self *= A
        self.obj *= A.obj
        self.char_Mv = False
        self.characterise_Mv()
        return(self)

    def __div_ab__(self,A):  # self /= A
        if isinstance(A,Mv):
            self *= A.inv()
        else:
            self *= S(1)/A
        return

    def __div__(self, A):
        if isinstance(A,Mv):
            return self * A.inv()
        else:
            return self * (S(1)/A)

    def __truediv__(self, A):
        if isinstance(A,Mv):
            return self * A.inv()
        else:
            return self * (S(1)/A)

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter
        return Printer().doprint(self)

    def __repr__(self):
        return str(self)

    def raw_str(self):
        return self.Mv_str(raw=True)

    def format_mv_coef(self,coef,base,first=False,latex_flg=False):

        coef_str = str(coef)
        base_str = str(base)
        add_flg  = isinstance(coef, Add)

        if base_str == '1':  # Case of scalar part of multivector
            if coef_str[0] == '-' or first:
                return coef_str
            else:
                return '+' + coef_str

        if coef_str == '1':
            if first:
                return base_str
            else:
                return '+' + base_str


        if coef_str == '-1':
            return '-' + base_str

        if latex_flg:
            if isinstance(coef, Add):
                if coef_str[0] == '-':
                    coef_str = r'-\left (' + str(-coef) + r'\right )'
                else:
                    coef_str = r'\left (' + coef_str + r'\right )'
        else:
            if isinstance(coef, Add):
                if coef_str[0] == '-':
                    coef_str = '-(' + str(-coef) + ')'
                else:
                    coef_str = '(' + coef_str + ')'

        if not first:
            if len(coef_str) > 0:
                if coef_str[0] != '-':
                    coef_str = '+' + coef_str
            else:
                coef_str = '+' + coef_str

        if latex_flg:
            return coef_str + base_str
        else:
            return coef_str + '*' + base_str


    def Mv_str(self,raw=False):
        if self.i_grade == 0:
            return str(self.obj)

        if self.is_blade_rep or self.Ga.is_ortho:
            grade_dict = self.Ga.blades_to_grades_dict
            key_sort = self.Ga.blade_key
        else:
            grade_dict = self.Ga.bases_to_grades_dict
            key_sort = self.Ga.base_key

        if not raw:
            self.obj = metric.Simp.apply(self.obj)

        (coefs, bases) = metric.linear_expand(self.obj)
        terms = list(zip(coefs,bases))
        terms_sorted = sorted(terms,key=key_sort)

        s = ''
        first = True

        if printer.GaPrinter.fmt == 1:  # Multivector printed on one line

            for (coef,base) in terms_sorted:

                s += self.format_mv_coef(coef,base,first=first)
                first = False

            return s

        elif printer.GaPrinter.fmt == 2:  # Multivector printed one grade per line

            grade_dict[1] = 1
            old_grade = -1

            for (coef,base) in terms_sorted:

                new_grade = grade_dict[base]

                if new_grade != old_grade and not first:
                    s += '\n'
                    old_grade = new_grade

                term_str = self.format_mv_coef(coef, base, first = first).strip()
                first = False
                s += term_str

            return s

        elif printer.GaPrinter.fmt == 3:  # Multivector printed one base per line

            for (coef,base) in terms_sorted:

                term_str = self.format_mv_coef(coef, base, first = first).strip()
                first = False
                s += term_str + '\n'

            return s[:-1]

    def raw_latex_str(self):
        return self.Mv_latex_str(raw=True)

    def Mv_latex_str(self,raw=False):
        if self.i_grade == 0:
            return str(self.obj)

        if self.is_blade_rep or self.Ga.is_ortho:
            grade_dict = self.Ga.blades_to_grades_dict
            key_sort = self.Ga.blade_key
        else:
            grade_dict = self.Ga.bases_to_grades_dict
            key_sort = self.Ga.base_key

        if not raw:
            self.obj = metric.Simp.apply(self.obj)

        (coefs, bases) = metric.linear_expand(self.obj)
        terms = list(zip(coefs,bases))
        terms_sorted = sorted(terms,key=key_sort)

        first = True
        s = ''

        if printer.GaLatexPrinter.fmt == 1:  # Multivector printed on one line

            for (coef,base) in terms_sorted:

                s += self.format_mv_coef(coef,base,first=first,latex_flg=True)
                first = False

            return s

        elif printer.GaLatexPrinter.fmt == 2:  # Multivector printed one grade per line

            grade_dict[1] = 1
            old_grade = -1

            for (coef,base) in terms_sorted:

                new_grade = grade_dict[base]

                if new_grade != old_grade and not first:
                    s += r' \\ '
                    old_grade = new_grade

                term_str = self.format_mv_coef(coef, base, first = first,latex_flg=True).strip()
                first = False
                s += term_str

            return r'\begin{array}{l} ' + s + r'\end{array} '

        elif printer.GaLatexPrinter.fmt == 3:  # Multivector printed one base per line

            for (coef,base) in terms_sorted[:-1]:

                term_str = self.format_mv_coef(coef, base, first = first,latex_flg=True).strip()
                first = False
                s += term_str + r' \\ '

            if len(terms_sorted) > 1:
                (coef,base) = terms_sorted[-1]
                term_str = self.format_mv_coef(coef, base, first = first,latex_flg=True).strip()
                s += term_str

            return r'\begin{array}{l} ' + s[:-4] + r'\end{array} '

    def __xor__(self, A):  # Mv wedge (^) product

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            return Mv(A * self.obj, ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In ^ operation Mv arguments are not from same geometric algebra')

        if isinstance(A, (Dop, Sdop, Pdop)):
            if isinstance(A, Pdop):
                A = A.promote().promote()
            if isinstance(A, Sdop):
                A = A.promote()
            return Dop.__rxor__(A, self)

        if self.is_scalar():
            return self * A

        self = self.blade_rep()
        A = A.blade_rep()
        self_W_A = self.Ga.wedge(self.obj, A.obj)
        self_W_A = Mv(self_W_A, ga=self.Ga)
        return self_W_A

    def __rxor__(self, A):  # Mv wedge (^) product
        if not isinstance(A, Mv):
            return Mv(A * self.obj, ga=self.Ga)
        else:
            return A * self

    def __or__(self, A):  # Mv dot (|) product
        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            return Mv(ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In | operation Mv arguments are not from same geometric algebra')

        self.Ga.dot_mode = '|'

        if isinstance(A, (Dop, Sdop, Pdop)):
            if isinstance(A, Pdop):
                A = A.promote().promote()
            if isinstance(A, Sdop):
                A = A.promote()
            return Dop.__ror__(A, self)

        self = self.blade_rep()
        if self.is_scalar() or A.is_scalar():
            return S(0)
        A = A.blade_rep()
        self_dot_A = Mv(self.Ga.dot(self.obj, A.obj), ga=self.Ga)
        return self_dot_A

    def __ror__(self, A):  # dot (|) product
        if not isinstance(A, Mv):
            return Mv(ga=self.Ga)
        else:
            return A | self

    def __pow__(self,n):  # Integer power operator
        if not isinstance(n,int):
            raise ValueError('!!!!Multivector power can only be to integer power!!!!')

        result = S(1)
        for x in range(n):
            result *= self
        return result

    def __lshift__(self, A): # anti-comutator (<<)
        return half * (self * A + A * self)

    def __rshift__(self, A): # comutator (>>)
        return half * (self * A - A * self)

    def __rlshift__(self, A): # anti-comutator (<<)
        return half * (A * self + self * A)

    def __rrshift__(self, A): # comutator (>>)
        return half * (A * self - self * A)

    def __lt__(self, A):  # left contraction (<)

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):  # sympy scalar
            return Mv(A * self.obj, ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In < operation Mv arguments are not from same geometric algebra')

        self.Ga.dot_mode = '<'

        if isinstance(A, Dop):
            return A.Mul(self, A, op='<')

        self = self.blade_rep()
        A = A.blade_rep()
        """
        if A.is_scalar():
            if self.is_scalar():
                return self.obj * A.obj
            else:
                return S(0)
        """

        self_lc_A = Mv(self.Ga.dot(self.obj, A.obj), ga=self.Ga)
        return self_lc_A

    def __gt__(self, A):  # right contraction (>)

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):  # sympy scalar
            return self.Ga.mv(A * self.scalar())

        if self.Ga.name != A.Ga.name:
            raise ValueError('In > operation Mv arguments are not from same geometric algebra')

        self.Ga.dot_mode = '>'

        if isinstance(A, Dop):
            return A.Mul(self, A, op='>')

        self = self.blade_rep()
        A = A.blade_rep()
        """
        if self.is_scalar():
            if A.is_scalar():
                return self.obj * A.obj
            else:
                return S(0)
        """

        self_rc_A = Mv(self.Ga.dot(self.obj, A.obj), ga=self.Ga)
        return self_rc_A

    def collect(self,deep=False):
        """
        # group coeffients of blades of multivector
        # so there is only one coefficient per grade
        self.obj = expand(self.obj)
        if self.is_blade_rep or Mv.Ga.is_ortho:
            c = self.Ga.blades_lst
        else:
            c = self.Ga.bases_lst
        self.obj = self.obj.collect(c)
        return self
        """
        coefs, bases = metric.linear_expand(self.obj)
        obj_dict = {}
        for (coef, base) in zip(coefs, bases):
            if base in list(obj_dict.keys()):
                obj_dict[base] += coef
            else:
                obj_dict[base] = coef
        obj = 0
        for base in list(obj_dict.keys()):
            if deep:
                obj += collect(obj_dict[base])*base
            else:
                obj += obj_dict[base]*base
        self.obj = obj
        return(self)


    def is_scalar(self):
        grades = self.Ga.grades(self.obj)
        if len(grades) == 1 and grades[0] == 0:
            return True
        else:
            return False

    def is_vector(self):
        grades = self.Ga.grades(self.obj)
        if len(grades) == 1 and grades[0] == 1:
            return True
        else:
            return False

    def is_blade(self):  # True is self is blade, otherwise False
        # sets self.blade_flg and returns value

        if self.blade_flg is not None:
            return self.blade_flg
        else:
            if self.is_versor():
                if self.i_grade is not None:
                    self.blade_flg = True
                else:
                    self.blade_flg = False
            else:
                self.blade_flg = False
            return self.blade_flg

    def is_base(self):
        (coefs, _bases) = metric.linear_expand(self.obj)
        if len(coefs) > 1:
            return False
        else:
            return coefs[0] == ONE

    def is_versor(self):  # Test for versor (geometric product of vectors)
        """
        This follows Leo Dorst's test for a versor.
        Leo Dorst, 'Geometric Algebra for Computer Science,' p.533
        Sets self.versor_flg and returns value
        """

        if self.versor_flg is not None:
            return self.versor_flg
        self.characterise_Mv()
        self.versor_flg = False
        self_rev = self.rev()
        # see if self*self.rev() is a scalar
        test = self*self_rev
        if not test.is_scalar():
            return self.versor_flg
        # see if self*x*self.rev() returns a vector for x an arbitrary vector
        test = self * self.Ga.XOX * self.rev()
        self.versor_flg = test.is_vector()
        return self.versor_flg

    def is_zero(self):
        self = self.simplify()
        if self.obj == 0:
            return True
        return False

    def scalar(self):
        # return scalar part of multivector
        # as sympy expression
        return self.Ga.scalar_part(self.obj)

    def get_grade(self, r):
        # return r-th grade of multivector as
        # a multivector
        return Mv(self.Ga.get_grade(self.obj, r), ga=self.Ga)

    def components(self):
        (coefs, bases) = metric.linear_expand(self.obj)
        bases_lst = self.Ga.blades_lst
        cb = list(zip(coefs, bases))
        cb = sorted(cb, key=lambda x: self.Ga.blades_lst0.index(x[1]))
        terms = []
        for (coef, base) in cb:
            terms.append(self.Ga.mv(coef * base))
        return terms

    def get_coefs(self, grade):
        (coefs, bases) = metric.linear_expand(self.obj)
        bases_lst = self.Ga.blades_lst
        cb = list(zip(coefs, bases))
        cb = sorted(cb, key=lambda x: self.Ga.blades[grade].index(x[1]))
        (coefs, bases) = zip(*cb)
        return list(coefs)

    def blade_coefs(self, blade_lst=None):
        """
        For a multivector, A, and a list of basis blades, blade_lst return
        a list (sympy expressions) of the coefficients of each basis blade
        in blade_lst
        """

        if blade_lst is None:
            blade_lst = [self.Ga.mv(ONE)] + self.Ga.mv_blades_lst

        #print 'Enter blade_coefs blade_lst =', blade_lst, type(blade_lst), [i.is_blade() for i in blade_lst]

        for blade in blade_lst:
            if not blade.is_base() or not blade.is_blade():
                raise ValueError("%s expression isn't a basis blade" % blade)
        blade_lst = [x.obj for x in blade_lst]
        (coefs, bases) = metric.linear_expand(self.obj)
        coef_lst = []
        for blade in blade_lst:
            if blade in bases:
                coef_lst.append(coefs[bases.index(blade)])
            else:
                coef_lst.append(ZERO)
        return coef_lst

    def proj(self, bases_lst):
        """
        Project multivector onto a given list of bases.  That is find the
        part of multivector with the same bases as in the bases_lst.
        """
        bases_lst = [x.obj for x in bases_lst]
        (coefs, bases) = metric.linear_expand(self.obj)
        obj = 0
        for (coef, base) in zip(coefs, bases):
            if base in bases_lst:
                obj += coef * base
        return Mv(obj, ga=self.Ga)

    def dual(self):
        mode = ga.Ga.dual_mode_value
        sign = S(1)
        if '-' in mode:
            sign = -sign
        if 'Iinv' in mode:
            I = self.Ga.i_inv
        else:
            I = self.Ga.i
        if mode[0] == '+' or mode[0] == '-':
            return sign * I * self
        else:
            return sign * self * I

    def even(self):
        # return even parts of multivector
        return Mv(self.Ga.even_odd(self.obj, True), ga=self.Ga)

    def odd(self):
        # return odd parts of multivector
        return Mv(self.Ga.even_odd(self.obj, False), ga=self.Ga)

    def rev(self):
        self = self.blade_rep()
        return Mv(self.Ga.reverse(self.obj), ga=self.Ga)

    def diff(self, coord):
        Dself = Mv(ga=self.Ga)
        if self.Ga.coords is None:
           Dself.obj = diff(self.obj, coord)
           return Dself
        elif coord not in self.Ga.coords:
            if self.Ga.par_coords is None:
                Dself.obj = diff(self.obj, coord)
            elif coords not in self.Ga.par_coords:
                Dself.obj = diff(self.obj, coord)
            else:
                Dself.obj = diff(self.obj, coord)
                for x_coord in self.Ga.coords:
                    f = self.Ga.par_coords[x_coord]
                    if f != S(0):
                        tmp1 = self.Ga.pDiff(self.obj, x_coord)
                        tmp2 = diff(f, coord)
                        Dself.obj += tmp1 * tmp2
            Dself.characterise_Mv()
            return Dself
        else:
            Dself.obj = self.Ga.pDiff(self.obj, coord)
            Dself.characterise_Mv()
            return Dself

    def pdiff(self, var):
        return Mv(self.Ga.pDiff(self.obj, var), ga=self.Ga)

    def Grad(self, coords, mode='*', left=True):
        """
        Returns various derivatives (*,^,|,<,>) of multivector functions
        with respect to arbitrary coordinates, 'coords'.  This would be
        used where you have a multivector function of both the basis
        coordinate set and and auxilliary coordinate set.  Consider for
        example a linear transformation in which the matrix coefficients
        depend upon the manifold coordinates, but the vector being
        transformed does not and you wish to take the divergence of the
        linear transformation with respect to the linear argument.
        """
        return Mv(self.Ga.Diff(self, mode, left, coords=coords), ga=self.Ga)

    def exp(self, hint='-'):  # Calculate exponential of multivector
        """
        Only works if square of multivector is a scalar.  If square is a
        number we can determine if square is > or < zero and hence if
        one should use trig or hyperbolic functions in expansion.  If
        square is not a number use 'hint' to determine which type of
        functions to use in expansion
        """
        self = self.blade_rep()
        self_sq = self * self
        if self_sq.is_scalar():
            sq = simplify(self_sq.obj)  # sympy expression for self**2
            if sq == S(0):  # sympy expression for self**2 = 0
                return self + S(1)
            (coefs,bases) = metric.linear_expand(self.obj)
            if len(coefs) == 1:  # Exponential of scalar * base
                base = bases[0]
                base_Mv = self.Ga.mv(base)
                base_sq = (base_Mv*base_Mv).scalar()
                if hint == '-': # base^2 < 0
                    base_n = sqrt(-base_sq)
                    return self.Ga.mv(cos(base_n*coefs[0]) + sin(base_n*coefs[0])*(bases[0]/base_n))
                else:  # base^2 > 0
                    base_n = sqrt(base_sq)
                    return self.Ga.mv(cosh(base_n*coefs[0]) + sinh(base_n*coefs[0])*(bases[0]/base_n))
            if sq.is_number:  # Square is number, can test for sign
                if sq > S(0):
                    norm = sqrt(sq)
                    value = self.obj / norm
                    tmp = Mv(cosh(norm) + sinh(norm) * value, ga=self.Ga)
                    tmp.is_blade_rep = True
                    return tmp
                else:
                    norm = sqrt(-sq)
                    value = self.obj / norm
                    tmp = Mv(cos(norm) + sin(norm) * value, ga=self.Ga)
                    tmp.is_blade_rep = True
                    return tmp
            else:
                if hint == '+':
                    norm = simplify(sqrt(sq))
                    value = self.obj / norm
                    tmp = Mv(cosh(norm) + sinh(norm) * value, ga=self.Ga)
                    tmp.is_blade_rep = True
                    return tmp
                else:
                    norm = simplify(sqrt(-sq))
                    value = self.obj / norm
                    obj = cos(norm) + sin(norm) * value
                    tmp = Mv(cos(norm) + sin(norm) * value, ga=self.Ga)
                    tmp.is_blade_rep = True
                    return tmp
        else:
            raise ValueError('"' + str(self) + '**2" is not a scalar in exp.')

    def set_coef(self, igrade, ibase, value):
        if self.blade_rep:
            base = self.Ga.blades[igrade][ibase]
        else:
            base = self.Ga.bases[igrade][ibase]
        (coefs, bases) = metric.linear_expand(self.obj)
        bases_lst = list(bases)  # python 2.5
        if base in bases:
            self.obj += (value - coefs[bases_lst.index(base)]) * base
        else:
            self.obj += value * base
        return

    def Fmt(self, fmt=1):
        """
        Set format for printing of multivectors -

            fmt = 1 - One multivector per line
            fmt = 2 - One grade per line
            fmt = 3 - one base per line

        Usage for multivector A example is -

            A.Fmt('2','A')

        output is

            'A = '+str(A)

        with one grade per line.  Works for both standard printing and
        for latex.
        """

        if printer.GaLatexPrinter.latex_flg:
            printer.GaLatexPrinter.fmt = fmt
            latex_str = printer.GaLatexPrinter.latex(self)
            printer.GaLatexPrinter.fmt = 1
            return latex_str
        else:
            printer.GaPrinter.fmt = fmt
            text_str = str(self)
            printer.GaPrinter.fmt = 1
            return text_str

        return

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        #if r'\begin{align*}' not in latex_str:
        if r'\begin{aligned}' not in latex_str:
            if self.title is None:
                latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
            else:
                latex_str = r'\begin{equation*} ' + self.title + ' = ' + latex_str + r' \end{equation*}'
        else:
            if self.title is not None:
                latex_str = latex_str.replace('&',' ' + self.title + ' =&',1)
        return latex_str

    def odot(self,dot_flg=True):
        new_self = copy.deepcopy(self)
        new_self.dot_flg = dot_flg
        return new_self

    def norm2(self):
        return self&self.rev()

    def norm(self, hint='+'):
        if hint == '+':
            return sqrt(self.norm2())
        else:
            return sqrt(Abs(self.norm2()))

    def inv(self):
        if self.is_scalar():  # self is a scalar
            return self.Ga.mv(S(1)/self.obj)
        self_sq = self * self
        if self_sq.is_scalar():  # self*self is a scalar
            """
            if self_sq.scalar() == S(0):
                raise ValueError('!!!!In multivector inverse, A*A is zero!!!!')
            """
            return (S(1)/self_sq.obj)*self
        self_rev = self.rev()
        self_self_rev = self * self_rev
        if(self_self_rev.is_scalar()): # self*self.rev() is a scalar
            """
            if self_self_rev.scalar() == S(0):
                raise ValueError('!!!!In multivector inverse A*A.rev() is zero!!!!')
            """
            return (S(1)/self_self_rev.obj) * self_rev
        raise TypeError('In inv() for self =' + str(self) + 'self, or self*self or self*self.rev() is not a scalar')

    def func(self, fct):  # Apply function, fct, to each coefficient of multivector
        (coefs, bases) = metric.linear_expand(self.obj)
        s = S(0)
        for (coef, base) in zip(coefs, bases):
            s += fct(coef) * base
        fct_self = Mv(s, ga=self.Ga)
        fct_self.characterise_Mv()
        return fct_self

    def trigsimp(self):
        return self.func(trigsimp)

    def simplify(self, modes=simplify):
        (coefs, bases) = metric.linear_expand(self.obj)
        obj = S(0)
        if isinstance(modes, list) or isinstance(modes, tuple):
            for (coef, base) in zip(coefs, bases):
                for mode in modes:
                    coef = mode(coef)
                obj += coef * base
        else:
            for (coef, base) in zip(coefs, bases):
                obj += modes(coef) * base
        self.obj = obj
        return self

    def subs(self, d):
        # For each scalar coef of the multivector apply substitution argument d
        (coefs, bases) = metric.linear_expand(self.obj)
        obj = S(0)
        for (coef, base) in zip(coefs, bases):
            obj += coef.subs(d) * base
        #self.obj = obj
        #return self
        return self.Ga.mv(obj)

    def expand(self):
        coefs,bases = metric.linear_expand(self.obj)
        new_coefs = []
        for coef in coefs:
            new_coefs.append(expand(coef))
        obj = 0
        for coef,base in zip(new_coefs,bases):
            obj += coef * base
        self.obj = obj
        return self

    def list(self):
        (coefs, bases) = metric.linear_expand(self.obj)
        indexes = []
        key_coefs = []
        for (coef, base) in zip(coefs, bases):
            if base in self.Ga.basis:
                index = self.Ga.basis.index(base)
                key_coefs.append((coef, index))
                indexes.append(index)

        for index in self.Ga.n_range:
            if index not in indexes:
                key_coefs.append((S(0), index))

        key_coefs = sorted(key_coefs, key=itemgetter(1))
        coefs = [x[0] for x in key_coefs]
        return coefs

    def grade(self, r=0):
        return self.get_grade(r)

    def pure_grade(self):
        """
        For pure grade return grade.  If not pure grade return negative
        of maximum grade
        """
        self.characterise_Mv()
        if self.i_grade is not None:
            return self.i_grade
        return -self.grades[-1]

def compare(A,B):
    """
    Determine is B = c*A where c is a scalar.  If true return c
    otherwise return 0.
    """
    if isinstance(A, Mv) and isinstance(B, Mv):
        Acoefs, Abases = metric.linear_expand(A.obj)
        Bcoefs, Bbases = metric.linear_expand(B.obj)
        if len(Acoefs) != len(Bcoefs):
            return 0
        if Abases != Bbases:
            return 0
        if Bcoefs[0] != 0 and Abases[0] == Bbases[0]:
            c = simplify(Acoefs[0]/Bcoefs[0])
            print('c =',c)
        else:
            return 0
        for acoef,abase,bcoef,bbase in zip(Acoefs[1:],Abases[1:],Bcoefs[1:],Bbases[1:]):
            print(acoef,'\n',abase,'\n',bcoef,'\n',bbase)
            if bcoef != 0 and abase == bbase:
                print('c-a/b =',simplify(c-(acoef/bcoef)))
                if simplify(acoef/bcoef) != c:
                    return 0
                else:
                    pass
            else:
                return 0
        return c
    else:
        raise TypeError('In compare both arguments are not multivectors\n')

#################### Partial Derivative Operator Class #################

class Pdop(object):
    """
    Partial derivative class for multivectors.  The partial derivatives
    are of the form

        \partial_{i_{1}...i_{n}} =
            \partial^{i_{1}+...+i_{n}}/\partial{x_{1}^{i_{1}}}...\partial{x_{n}^{i_{n}}}.

    If i_{j} = 0 then the partial derivative does not contain the x^{i_{j}}
    coordinate.

    The partial derivative is represented by a dictionary with coordinates
    for keys and key value are the number of times one differentiates with
    respect to the key.
    """

    ga = None

    init_slots = {'ga': (None, 'Associated geometric algebra')}

    @staticmethod
    def setGa(ga):
        Pdop.Ga = ga
        return

    @staticmethod
    def compare_pdiffs(Pdop1, Pdop2):  # compare two Pdops

        if isinstance(Pdop1, tuple):
            pdop1 = Pdop1[1]
        else:
            pdop1 = Pdop1

        if isinstance(Pdop2, tuple):
            pdop2 = Pdop2[1]
        else:
            pdop2 = Pdop2


        if pdop1.order > pdop2.order:
            return 1
        if pdop1.order < pdop2.order:
            return -1

        keys1 = list(pdop1.pdiffs.keys())
        keys2 = list(pdop2.pdiffs.keys())
        lkeys1 = len(keys1)
        lkeys2 = len(keys2)

        if lkeys1 == lkeys2:
            s1 = ''.join([str(pdop1.Ga.coords.index(x)) for x in keys1])
            s2 = ''.join([str(pdop1.Ga.coords.index(x)) for x in keys2])
            if s1 < s2:
                return -1
            else:
                return 1
        else:
            if lkeys1 < lkeys2:
                return 1
            else:
                return -1

        return 0

    def __eq__(self,A):
        if isinstance(A, Pdop) and self.Ga.name == A.Ga.name and self.pdiffs == A.pdiffs:
            return True
        else:
            if len(self.pdiffs) == 0 and A == S(1):
                return True
            return False

    def __init__(self, *kargs, **kwargs):
        """
        The partial differential operator is a partial derivative with
        respect to a set of real symbols (variables).  The allowed
        variables are in two lists.  self.Ga.coords is a list of the
        coordinates associated with the geometric algebra.  self.Ga.auxvars
        is a list of auxiallary symbols that have be added to the geometric
        algebra using the member function Ga.AddVars(self,auxvars).

        The data structure of a Pdop is the dictionary self.pdiffs where
        the keys are the variables of differentiation and the values are the
        order of differentiation of with respect to the variable.

        Examples:
            First dervative with respect to x -
                PDx = Pdop(x,ga=o3d)
            Second dervative with respect to x -
                P2Dx = Pdop({x:2},ga=o3d)
            First derivative with respect to both x and y -
                P2Dx = Pdop({x:1,y:1},ga=o3d)
         """

        kwargs = metric.test_init_slots(Pdop.init_slots, **kwargs)

        self.Ga = kwargs['ga']  # Associated geometric algebra
        self.order = 0

        if self.Ga is None:
            raise ValueError('In Pdop.__init__ self.Ga must be defined.')

        if kargs[0] is None:  # Pdop is the identity (1)
            self.pdiffs = {}
        elif isinstance(kargs[0], dict):  # Pdop defined by dictionary
            self.pdiffs = kargs[0]
        elif isinstance(kargs[0],Symbol):  # First order derivative with respect to symbol
            self.pdiffs = {kargs[0]:1}
        else:
            raise ValueError('In pdop kargs = ', str(kargs))

        for x in list(self.pdiffs.keys()):  # self.order is total number of differentiations
            self.order += self.pdiffs[x]

    def factor(self):
        """
        If partial derivative operator self.order > 1 factor out first
        order differential operator.  Needed for application of partial
        derivative operator to product of sympy expression and partial
        differential operator.  For example if D = Pdop({x:3}) then

            (Pdop({x:2}),Pdop({x:1})) = D.factor()
        """
        if self.order == 1:
            return S(0), self
        else:
            x = list(self.pdiffs.keys())[0]
            self.order -= 1
            n = self.pdiffs[x]
            if n == 1:
                del self.pdiffs[x]
            else:
                self.pdiffs[x] -= 1
            return self, self.Ga.Pdiffs[x]

    def __call__(self, arg):
        """
        Calculate nth order partial derivative (order defined by
        self) of Mv, Dop, Sdop or sympy expression
        """
        if self.pdiffs == {}:
            return arg  # result is Pdop identity (1)

        if isinstance(arg, Pdop):  # arg is Pdop
            if self.Ga.name != arg.Ga.name:
                raise ValueError('In Pdop.__call__ arguments do not belong to same geometric algebra.')
            elif arg.pdiffs == {}:  # arg is one
                #return self
                return S(0)  # derivative is zero
            else:  # arg is partial derivative
                pdiffs = copy.copy(arg.pdiffs)
                for key in self.pdiffs:
                    if key in pdiffs:
                        pdiffs[key] += self.pdiffs[key]
                    else:
                        pdiffs[key] = self.pdiffs[key]
            return Pdop(pdiffs,ga=self.Ga)  # result is Pdop

        elif isinstance(arg, Dop):  # arg is Dop
            terms = []
            for x in arg.terms:
                px0 = self(x[0])
                if px0 != 0:
                    terms.append((self(x[0]),x[1]))
                terms.append((x[0],self(x[1])))
            return Dop(terms,ga=arg.Ga)

        elif isinstance(arg, Mv):  # arg is multivector
            for x in self.pdiffs:
                for i in range(self.pdiffs[x]):
                    arg = self.Ga.pDiff(arg, x)
            return arg  # result is multivector

        elif isinstance(arg, lt.Lt):
            SD_lt = {}
            for key in arg.lt_dict:  # Apply SDop to linear transformation dictionary
                SD_lt[key] = self(arg.lt_dict[key])
            return lt.Lt(SD_lt,ga=arg.Ga)

        elif isinstance(arg, (Expr, Symbol, numbers.Number)):  # arg is sympy expression
            for x in self.pdiffs:
                arg = diff(arg,x,self.pdiffs[x])
            return arg  # derivative is sympy expression

        elif isinstance(arg, list):  # arg is list of tuples (coef, partial derivative)
            D = copy.deepcopy(self)
            terms = copy.deepcopy(arg)
            while True:
                D, D0 = D.factor()
                k = 0
                for term in terms:
                    dc = D0(term[0])
                    pd = D0(term[1])
                    #print 'D0, term, dc, pd =', D0, term, dc, pd
                    tmp = []
                    if dc != 0:
                        tmp.append((dc,term[1]))
                    if pd != 0 :
                        tmp.append((term[0],pd))
                    terms[k] = tmp
                    k += 1
                terms = [i for o in terms for i in o]  # flatten list one level
                if D == 0:
                    break
            terms = Sdop.consolidate_coefs(terms)
            return terms  # result is list of tuples (coef, partial derivative)
        elif isinstance(arg, Sdop):  # arg is scalar differential operator
            if self.Ga != arg.Ga:
                raise ValueError('In Pdop.__call__ self.Ga != arg.Ga.')
            return self(arg.terms)  # result is list of tuples (coef, partial derivative)
        else:
            raise ValueError('In Pdop.__call__ type(arg) = ' + str(type(arg)) + ' not allowed.')

    def __mul__(self, arg):  # functional product of self and arg (self*arg)
        if isinstance(arg, (Sdop, Dop)):
            if isinstance(arg, Sdop):
                self = self.promote()
            self = self.promote()
            return self*arg
        return self(arg)

    def __rmul__(self, arg):
        if isinstance(arg, Mv):
            return arg*self.promote().promote()
        return arg*self.promote()

    def __add__(self, arg):
        if isinstance(arg, Pdop):
            return self.promote()+arg.promote()
        elif isinstance(arg, Sdop):
            return self.promote()+arg
        elif isinstance(arg, (Dop, Mv)):
            if isinstance(arg, Mv) and arg.is_scalar():
                return self.promote()+arg.obj
            return self.promote().promote()+arg
        else:
            return self.promote()+arg

    def __sub__(self,arg):
        if isinstance(arg, Pdop):
            return self.promote()-arg.promote()
        elif isinstance(arg, Sdop):
            return self.promote()-arg
        elif isinstance(arg, (Dop, Mv)):
            if isinstance(arg, Mv) and arg.is_scalar():
                return self.promote()-arg.obj
            return self.promote().promote()-arg
        else:
            return self.promote()-arg

    def __radd__(self, arg):
        return self.promote()+arg

    def __rsub__(self, arg):
        return self.promote()-arg

    def __neg__(self):
        return -self.promote()

    def __pos__(self):
        return self.promote()

    def Pdop_str(self):
        if self.order == 0:
            return '1'
        s = 'D'
        for x in self.pdiffs:
            s += '{' + str(x) + '}'
            n = self.pdiffs[x]
            if n > 1:
                s += '^' + str(n)
        return s

    def promote(self):
        return Sdop([1],[self],ga=self.Ga)

    def Pdop_latex_str(self):
        if self.order == 0:
            return '1'
        s = r'{\displaystyle\frac{\partial'
        if self.order > 1:
            s += '^{' + str(self.order) + '}'
        s += '}{'
        keys = list(self.pdiffs.keys())
        keys.sort(key=(self.Ga.coords + keys).index)
        for key in keys:
            i = self.pdiffs[key]
            s += r'\partial ' + str(key)
            if i > 1:
                s += '^{' + str(i) + '}'
        s += '}}'
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



################ Scalar Partial Differential Operator Class ############

class Sdop(object):
    """
    Scalar differential operator is of the form (Einstein summation)

        D = c_{i}*D_{i}

    where the c_{i}'s are scalar coefficient (they could be functions)
    and the D_{i}'s are partial differential operators.
    """

    init_slots = {'ga': (None, 'Associated geometric algebra')}

    ga = None
    str_mode = False

    @staticmethod
    def setGa(ga):
        Sdop.Ga = ga
        Pdop.setGa(ga)
        return

    def TSimplify(self):
        new_terms = []
        for (coef, pdiff) in self.terms:
            new_terms.append((Simp.apply(coef), pdiff))
        self.terms = new_terms
        return

    def is_zero(self):
        for term in self.terms:
            if term[0] != 0:
                return False
        return True

    @staticmethod
    def consolidate_coefs(sdop):
        """
        Remove zero coefs and consolidate coefs with repeated pdiffs.
        """
        if isinstance(sdop, Sdop):
            terms = sdop.terms
        else:
            terms = sdop

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
        new_terms = list(zip(new_coefs, new_pdiffs))

        if isinstance(sdop, Sdop):
            return Sdop(new_terms, ga=sdop.Ga)
        else:
            return new_terms

    def simplify(self, modes=simplify):
        coefs, pdiffs = zip(*self.terms)

        coefs = list(coefs)
        pdiffs = list(pdiffs)
        new_coefs = []
        for coef in coefs:
            new_coefs.append(metric.apply_function_list(modes,coef))
        self.terms = list(zip(new_coefs,pdiffs))
        return self

    def sort_terms(self):

        if len(self.terms) > 1:
            self.terms = sorted(self.terms,key=functools.cmp_to_key(Pdop.compare_pdiffs))

        return

    def Sdop_str(self):

        self.sort_terms()

        s = ''
        for (coef, pdop) in self.terms:
            pd_str = str(pdop)

            if coef == S(1):
                s += pd_str
            elif coef == S(-1):
                s += '-' + pd_str
            else:
                if isinstance(coef, Add):
                    s += '(' + str(coef) + ')*' + pd_str
                else:
                    s += str(coef) + '*' + pd_str
            s += ' + '

        s = s.replace('+ -','- ')
        s = s[:-3]
        if Sdop.str_mode:
            if len(self.terms) > 1 or isinstance(self.terms[0][0], Add):
                s = '(' + s + ')'
        return s

    def Sdop_latex_str(self):
        if len(self.terms) == 0:
            return '0'

        self.sort_terms()
        first = True

        s = ''
        for (coef, pdop) in self.terms:

            pdop_str = str(pdop)

            coef_str = printer.format_coef(coef,first=first,latex_flg=True)
            first = False

            if pdop_str == '1':
                if coef_str[-1] == ')' and coef_str[0] != '-':
                    s += coef_str[1:-1]
                else:
                    s += coef_str
            else:
                s += coef_str + pdop_str

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

    def __init__(self, *kargs, **kwargs):
        """
        The scalar differential operator structure is of the form
        (Einstein summation)

            D = c_{i}D_{i}

        where the c_{i}'s are scalar coefficients and the D_{i}'s are
        partial differential operator (class Pdop).  D is stored in
        the structure self.terms = [(c_{1},D_{1}),(c_{2},D_{2}),...].
        """
        kwargs = metric.test_init_slots(Sdop.init_slots, **kwargs)

        self.Ga = kwargs['ga']  # Associated geometric algebra (coords)

        if self.Ga is None:
            raise ValueError('In Sdop.__init__ self.Ga must be defined.')

        if len(kargs[0]) == 0:  # identity Dop
            self.terms = [(S(1), self.Ga.Pdop_identity)]
        elif len(kargs[0]) == 1 and isinstance(kargs[0],Symbol):  # Simple Pdop of order 1
            self.terms = [(S(1), self.Ga.pdop(kargs[0]))]
        else:
            if len(kargs) == 2 and isinstance(kargs[0],list) and isinstance(kargs[1],list):
                if len(kargs[0]) != len(kargs[1]):
                    raise ValueError('In Sdop.__init__ coefficent list and Pdop list must be same length.')
                self.terms = list(zip(kargs[0],kargs[1]))
            elif len(kargs) == 1 and isinstance(kargs[0],list):
                self.terms = kargs[0]
            else:
                raise ValueError('In Sdop.__init__ length of kargs must be 1 or 2 kargs = '+str(kargs))

    def __call__(self, arg):
        if isinstance(arg, Sdop):
            if self.Ga != arg.Ga:
                raise ValueError('In Sdop.__call__  self.Ga != arg.Ga.')
            terms = []
            for (coef, pdiff) in self.terms:
                new_terms = pdiff(arg.terms)
                new_terms = [ (coef * x[0], x[1]) for x in new_terms]
                terms += new_terms
            return Sdop(terms, ga=self.Ga)
        elif isinstance(arg,lt.Lt):
            SD_lt = {}
            for key in arg.lt_dict:  # Apply SDop to linear transformation dictionary
                SD_lt[key] = self(arg.lt_dict[key])
            return lt.Lt(SD_lt,ga=arg.Ga)
        else:
            return sum([x[0] * x[1](arg) for x in self.terms])

    def __neg__(self):
        return Sdop([(-x[0], x[1]) for x in self.terms], ga=self.Ga)

    @staticmethod
    def Add(sdop1, sdop2):

        if isinstance(sdop1, Sdop) and isinstance(sdop2, Sdop):
            if sdop1.Ga != sdop2.Ga:
                raise ValueError('In Sdop.Add sdop1.Ga != sdop2.Ga.')
            coefs1, pdiffs1 = zip(*sdop1.terms)
            coefs2, pdiffs2 = zip(*sdop2.terms)

            coefs1 = list(coefs1)
            coefs2 = list(coefs2)

            pdiffs1 = list(pdiffs1)
            pdiffs2 = list(pdiffs2)

            pdiffs = pdiffs1 + [x for x in pdiffs2 if x not in pdiffs1]
            coefs = len(pdiffs) * [S(0)]

            for pdiff in pdiffs1:
                index = pdiffs.index(pdiff)
                coef = coefs1[pdiffs1.index(pdiff)]
                coefs[index] += coef

            for pdiff in pdiffs2:
                index = pdiffs.index(pdiff)
                coef = coefs2[pdiffs2.index(pdiff)]
                coefs[index] += coef

            sdop_sum = Sdop(coefs, pdiffs, ga=sdop1.Ga)

            if sdop_sum.is_zero():
                return 0

        elif isinstance(sdop1, Sdop):  # sdop1 is Sdop
            if isinstance(sdop2, Mv):
                return sdop1.promote()+sdop2
            else:
                coefs, pdiffs = zip(*sdop1.terms)
                coefs = list(coefs)
                pdiffs = list(pdiffs)
                if sdop1.Ga.Pdop_identity in pdiffs:
                    index = pdiffs.index(sdop1.Ga.Pdop_identity)
                    coefs[index] += sdop2
                else:
                    coefs.append(sdop2)
                    pdiffs.append(sdop1.Ga.Pdop_identity)

                sdop_sum = Sdop(coefs, pdiffs, ga=sdop1.Ga)

                if sdop_sum.is_zero():
                    return 0

        else:  # sdop2 is Sdop
            if isinstance(sdop1, Mv):
                if sdop1.is_scalar():
                    return Sdop.Add(sdop1.obj,sdop2)
                return Dop.Add(sdop1,sdop2.promote())
            else:
                coefs, pdiffs = zip(*sdop2.terms)
                coefs = list(coefs)
                pdiffs = list(pdiffs)
                if sdop2.Ga.Pdop_identity in pdiffs:
                    index = pdiffs.index(sdop2.Ga.Pdop_identity)
                    coefs[index] += sdop1
                else:
                    coefs.append(sdop1)
                    pdiffs.append(sdop2.Ga.Pdop_identity)

                sdop_sum = Sdop(coefs, pdiffs, ga=sdop2.Ga)

                if sdop_sum.is_zero():
                    return 0

        sdop_sum = Sdop.consolidate_coefs(sdop_sum)

        if sdop_sum.is_zero():
            return 0

        return sdop_sum

    def __eq__(self, sdop):
        if isinstance(sdop, Sdop):
            if self.Ga != sdop.Ga:
                return False
            self = Sdop.consolidate_coefs(self)
            sdop = Sdop.consolidate_coefs(sdop)
            if len(self.terms) != len(sdop.terms):
                return False
            if set(self.terms) != set(sdop.terms):
                return False
            return True
        else:
            return False

    def __add__(self, sdop):
        if isinstance(sdop,Pdop):
            return Sdop.Add(self,sdop.promote())
        return Sdop.Add(self, sdop)

    def __radd__(self, sdop):
        return Sdop.Add(sdop, self)

    def __add_ab__(self, sdop):
        if isinstance(sdop, Sdop):
            if self.Ga != sdop.Ga:
                raise ValueError('In Sdop.__add_ab__ self.Ga != sdop.Ga.')

            coefs, pdiffs = zip(*self.terms)
            pdiffs = list(pdiffs)
            coefs = list(coefs)

            for (coef, pdiff) in sdop.terms:
                if pdiff in pdiffs:
                    index = pdiffs.index(pdiff)
                    coefs[index] += coef
                else:
                    pdiffs.append(pdiff)
                    coefs.append(coef)
            self.term = list(zip(coefs, pdiffs))
            self = Sdop.consolidate_coefs(self)
            return

        elif isinstance(sdop, tuple):
            self.term.append(sdop)
            self = Dfop.consolidate_coefs(self)
            return

        else:
            self.terms.append((sdop, self.Ga.Pdop_identity))
            self = Sdop.consolidate_coefs(self)
            return

    def __sub__(self, sdop):
        sdop_sub = Sdop.Add(self, -sdop)
        return sdop_sub

    def __rsub__(self, sdop):
        sdop_sub = Sdop.Add(-self, sdop)
        if sdop_sub.is_zero():
            return 0
        return sdop_sub

    def __mul__(sdopl, sdopr):
        if isinstance(sdopr, (Sdop, Pdop, Dop)):
            if sdopl.Ga != sdopr.Ga:
                raise ValueError('In Sdop.__mul__ Sdop arguments are not from same geometric algebra')
        terms = []
        if isinstance(sdopr, Sdop):
            for (coef, pdiff) in sdopl.terms:
                Dsdopl = pdiff(sdopr.terms)  # list of terms
                Dsdopl = [(coef * x[0], x[1]) for x in Dsdopl]
                terms += Dsdopl
            product = Sdop(terms, ga=sdopl.Ga)
            return Sdop.consolidate_coefs(product)
        elif isinstance(sdopr, Pdop):
            return sdopl*sdopr.promote()
        elif isinstance(sdopr, Dop):
            return sdopl.promote()*sdopr
        else:
            return sdopl(sdopr)

    def __rmul__(self,sdop):
        terms = [(sdop * x[0], x[1]) for x in self.terms]
        return Sdop(terms, ga=self.Ga)

    def promote(self):
        coefs, pdiffs = zip(*self.terms)

        coefs = list(coefs)
        new_coefs = [ x*self.Ga.one for x in coefs]
        pdiffs = list(pdiffs)

        return Dop(new_coefs, pdiffs,ga=self.Ga)


################# Multivector Differential Operator Class ##############

class Dop(object):
    """
    Differential operator class for multivectors.  The operators are of
    the form

        D = D^{i_{1}...i_{n}}\partial_{i_{1}...i_{n}}

    where the D^{i_{1}...i_{n}} are multivector functions of the coordinates
    x_{1},...,x_{n} and \partial_{i_{1}...i_{n}} are partial derivative
    operators

        \partial_{i_{1}...i_{n}} =
             \partial^{i_{1}+...+i_{n}}/\partial{x_{1}^{i_{1}}}...\partial{x_{n}^{i_{n}}}.

    If * is any multivector multiplicative operation then the operator D
    operates on the multivector function F by the following definitions

        D*F = D^{i_{1}...i_{n}}*\partial_{i_{1}...i_{n}}F

    returns a multivector and

        F*D = F*D^{i_{1}...i_{n}}\partial_{i_{1}...i_{n}}

    returns a differential operator.  If the 'cmpflg' in the operator is
    set to 'True' the operation returns

        F*D = (\partial_{i_{1}...i_{n}}F)*D^{i_{1}...i_{n}}

    a multivector function.  For example the representation of the grad
    operator in 3d would be:

        D^{i_{1}...i_{n}} = [e__x,e__y,e__z]
        \partial_{i_{1}...i_{n}} = [(1,0,0),(0,1,0),(0,0,1)].

    See LaTeX documentation for definitions of operator algebraic
    operations +, -, *, ^, |, <, and >.
    """

    init_slots = {'ga': (None, 'Associated geometric algebra'),
                  'cmpflg': (False, 'Complement flag for Dop'),
                  'debug': (False, 'True to print out debugging information'),
                  'fmt_dop': (1, '1 for normal dop partial derivative formating')}

    ga = None

    @staticmethod
    def setGa(ga):  # set geometric algebra globally for all Dop's
        Dop.Ga = ga
        Sdop.setGa(ga)
        return

    @staticmethod
    def compare_terms(t1,t2):
        return Pdop.compare_pdiffs(t1[1],t2[1])

    @staticmethod
    def flatten_one_level(lst):
        return [inner for outer in lst for inner in outer]

    def __init__(self, *kargs, **kwargs):

        kwargs = metric.test_init_slots(Dop.init_slots, **kwargs)

        self.cmpflg = kwargs['cmpflg']  # Complement flag (default False)
        self.Ga = kwargs['ga']  # Associated geometric algebra
        self.dot_flg = False  # Flag for over dot notation

        if self.Ga is None:
            raise ValueError('In Dop.__init__ self.Ga must be defined.')

        self.dop_fmt = kwargs['fmt_dop']  # Partial derivative output format (default 1)
        self.title = None

        if len(kargs[0]) == 0:  # identity Dop
            self.terms = [(S(1),self.Ga.Pdop_identity)]
        else:
            if len(kargs) == 2:
                if len(kargs[0]) != len(kargs[1]):
                    raise ValueError('In Dop.__init__ coefficent list and Pdop list must be same length.')
                self.terms = list(zip(kargs[0],kargs[1]))
            elif len(kargs) == 1:
                if isinstance(kargs[0][0][0], Mv):  # Mv expansion [(Mv, Pdop)]
                    self.terms = kargs[0]
                elif isinstance(kargs[0][0][0], Sdop):  # Sdop expansion [(Sdop, Mv)]
                    coefs = []
                    pdiffs = []
                    for (sdop, mv) in kargs[0]:
                        for (coef, pdiff) in sdop.terms:
                            if pdiff in pdiffs:
                                index = pdiffs.index(pdiff)
                                coefs[index] += coef * mv
                            else:
                                pdiffs.append(pdiff)
                                coefs.append(coef * mv)
                    self.terms = list(zip(coefs, pdiffs))
                else:
                    raise ValueError('In Dop.__init__ kargs[0] form not allowed. kargs = ' + str(kargs))
            else:
                raise ValueError('In Dop.__init__ length of kargs must be 1 or 2.')

    def odot(self,dot_flg=True):
        new_self = copy.deepcopy(self)
        new_self.dot_flg = dot_flg
        return new_self

    def rop(self):  # Dop operates on right operands
        new_self = copy.deepcopy(self)
        new_self.cmpflg = False
        return new_self

    def lop(self):  # Dop operates on left operands
        new_self = copy.deepcopy(self)
        new_self.cmpflg = True
        return new_self

    def simplify(self, modes=simplify):
        """
        Simplify each multivector coefficient of a partial derivative
        """
        new_coefs = []
        new_pd = []
        for (coef, pd) in self.terms:
            tmp = coef.simplify(modes=modes)
            new_coefs.append(tmp)
            new_pd.append(pd)
        self.terms = list(zip(new_coefs, new_pd))
        return Dop(new_coefs, new_pd, ga=self.Ga, cmpflg=self.cmpflg)

    def consolidate_coefs(self):
        """
        Remove zero coefs and consolidate coefs with repeated pdiffs.
        """
        new_coefs = []
        new_pdiffs = []
        for (coef, pd) in self.terms:
            if isinstance(coef, Mv) and coef.is_scalar():
                coef = coef.obj
            if coef != S(0):
                if pd in new_pdiffs:
                    index = new_pdiffs.index(pd)
                    new_coefs[index] += coef
                else:
                    new_coefs.append(coef)
                    new_pdiffs.append(pd)

        self.terms = list(zip(new_coefs, new_pdiffs))
        return Dop(new_coefs, new_pdiffs, ga=self.Ga, cmpflg=self.cmpflg)


    def blade_rep(self):
        N = len(self.blades)
        coefs = N * [[]]
        bases = N * [0]
        for term in self.terms:
            for (coef, base) in metric.linear_expand(self.terms[0].obj, mode=False):
                index = self.blades.index(base)
                coefs[index] = coef
                bases[index] = base

    @staticmethod
    def Add(dop1, dop2):

        if isinstance(dop1, Dop) and isinstance(dop2, Dop):
            if dop1.Ga.name != dop2.Ga.name:
                raise ValueError('In Dop.Add Dop arguments are not from same geometric algebra')

            if dop1.cmpflg != dop2.cmpflg:
                raise ValueError('In Dop.Add complement flags have different values.')

            coefs1, pdiffs1 = zip(*dop1.terms)
            coefs2, pdiffs2 = zip(*dop2.terms)

            coefs1 = list(coefs1)
            coefs2 = list(coefs2)

            pdiffs1 = list(pdiffs1)
            pdiffs2 = list(pdiffs2)

            pdiffs = pdiffs1 + [x for x in pdiffs2 if x not in pdiffs1]
            coefs = len(pdiffs) * [S(0)]

            for pdiff in pdiffs1:
                index = pdiffs.index(pdiff)
                coef = coefs1[pdiffs1.index(pdiff)]
                coefs[index] += coef

            for pdiff in pdiffs2:
                index = pdiffs.index(pdiff)
                coef = coefs2[pdiffs2.index(pdiff)]
                coefs[index] += coef

            sum_dop = Dop(coefs, pdiffs, cmpflg=dop1.cmpflg, ga=dop1.Ga)

            if sum_dop.is_zero():
                return 0

            return sum_dop

        elif isinstance(dop1, Dop):  # dop1 is Dop:
            if isinstance(dop2, Pdop):
                return Dop.Add(dop1,dop2.promote().promote())
            elif isinstance(dop2, Sdop):
                return Dop.Add(dop1,dop2.promote())
            else:
                if not isinstance(dop2, Mv):
                    dop2 = dop1.Ga.mv(dop2)
                dop2 = Dop([dop2], [dop1.Ga.Pdop_identity], cmpflg=dop1.cmpflg, ga=dop1.Ga)
        elif isinstance(dop2, Dop):
            if isinstance(dop1, Pdop):
                return Dop.Add(dop1.promote().promote(),dop2)
            elif isinstance(dop1, Sdop):
                return Dop.Add(dop1.promote(),dop2)
            else:
                if not isinstance(dop1, Mv):  # dop2 is Dop
                    dop1 = dop2.Ga.mv(dop1)
                dop1 = Dop([dop1], [dop2.Ga.Pdop_identity], cmpflg=dop2.cmpflg, ga=dop2.Ga)
        return Dop.Add(dop1, dop2)

    def __add__(self, dop):
        return Dop.Add(self, dop)

    def __radd__(self, dop):
        return Dop.Add(dop, self)

    def __neg__(self):

        coefs, pdiffs = zip(*self.terms)

        coefs = list(coefs)
        pdiffs = list(pdiffs)

        coefs = [-x for x in coefs]

        neg = Dop(coefs, pdiffs, ga=self.Ga,
                  cmpflg=self.cmpflg)

        return neg

    def __sub__(self, dop):
        return Dop.Add(self, -dop)

    def __rsub__(self, dop):
        return Dop.Add(dop, -self)

    @staticmethod
    def Mul(dopl, dopr, op='*',diff=True):  # General multiplication of Dop's

        # Dop cmpflg is True if the Dop operates on the left argument and
        # False if the Dop operates on the right argument.

        # If diff = True and one agrument is not a Dop then operate with
        # partial derivatives in Dop argument on non Dop argument
        # otherwise just multiply the Dop coefficients by the non Dop
        # argument.  This is needed to implement the over dot methods.

        if isinstance(dopr, Pdop):
                dopr = dopr.promote()
        if isinstance(dopr, Sdop):
                dopr = dopr.promote()
        if isinstance(dopl, Pdop):
                dopl = dopl.promote()
        if isinstance(dopl, Sdop):
                dopl = dopl.promote()

        if isinstance(dopl, Dop) and isinstance(dopr, Dop):  # Both dopl and dopr are operators

            if dopl.Ga != dopr.Ga:
                raise ValueError('In Dop.Mul Dop arguments are not from same geometric algebra')
            if dopl.cmpflg != dopr.cmpflg:  # Cannot multiply left and right operators
                raise ValueError('In Dop.Mul Dop arguments do not have same cmplfg')

            if not dopl.cmpflg:  # dopl operates on right argument
                terms = []
                for (coef, pdiff) in dopl.terms:  # Apply each dopl term to dopr
                    Ddopl = pdiff(dopr.terms)  # list of terms
                    Ddopl = [(Mv.Mul(coef, x[0], op=op), x[1]) for x in Ddopl]
                    terms += Ddopl
                product = Dop(terms, ga=dopl.Ga)
            else:  # dopr operates on left argument
                terms = []
                for (coef, pdiff) in dopr.terms:  # Apply each dopr term to dopl
                    Ddopr = pdiff(dopl.terms)  # list of terms
                    Ddopr = [(Mv.Mul(x[0], coef, op=op), x[1]) for x in Ddopr]
                    terms += Ddopr
                product = Dop(terms, ga=dopr.Ga, cmpflg=True)

        else:  # Either dopl or dopr are differential operators (not both)

            if not isinstance(dopl, Dop):  # dopl is a scalar or Mv and dopr is Dop

                if isinstance(dopl, Mv) and dopl.Ga != dopr.Ga:
                    raise ValueError('In Dop.Mul Dop arguments are not from same geometric algebra')

                if dopr.cmpflg:  # dopr is a left operator
                    if not diff:  # Only multipy dopr coeficients by dopl
                        terms = [(Mv.Mul(dopl, x[0], op=op), x[1]) for x in dopr.terms]
                        return Dop(terms, ga=dopr.Ga, cmpflg=True)  # returns left Dop
                    else:  # Operate on dopl with dopr
                        product = sum([Mv.Mul(x[1](dopl), x[0], op=op) for x in dopr.terms])
                        return product
                else:  # dopr is a right operator only multiply dopr coefficients by dopl
                    terms = [(Mv.Mul(dopl, x[0], op=op), x[1]) for x in dopr.terms]
                    return Dop(terms, ga=dopr.Ga)

            else:  # dopr is a scalar or a Mv and dopl is a Dop

                if isinstance(dopr, Mv) and dopl.Ga != dopr.Ga:
                    raise ValueError('In Dop.Mul Dop arguments are not from same geometric algebra')

                if dopl.cmpflg:  # dopl is left operator
                    terms = [(Mv.Mul(x[0], dopr, op=op), x[1]) for x in dopl.terms]
                    product = Dop(terms, ga=dopl.Ga, cmpflg=True)  # dopl is a left operator only multiply dopl coefficients by dopr
                else:
                    if not diff:  # Only multiply dopl coefficients by dopr
                        terms = [(Mv.Mul(x[0], dopr, op=op), x[1]) for x in dopl.terms]
                        return Dop(terms, ga=dopr.Ga)  # returns Dop
                    else: # Operate on dopr with dopl
                        return sum([Mv.Mul(x[0], x[1](dopr), op=op) for x in dopl.terms])  # returns multivector

        if isinstance(product, Dop):
            product.consolidate_coefs()

        return product

    def TSimplify(self):
        new_terms = []
        for (coef, pdiff) in self.terms:
            new_terms.append((metric.Simp.apply(coef), pdiff))
        self.terms = new_terms
        return

    def __div__(self, a):
        if isinstance(a, (Mv, Dop)):
            raise TypeError('!!!!Can only divide Dop by sympy scalar expression!!!!')
        else:
            return (1/a) * self

    def __mul__(self, dopr):  # * geometric product

        """
        Effect of overdot on Dop __mul__

        Dop         *  Dop          = Dop
        Dop         *  Mv           = Mv
        Dop/Dot     *  Mv           = Dop/Dot
        Dop/Dot     *  Mv/Dot       = Mv
        """

        if not self.dot_flg:
            return Dop.Mul(self, dopr, op='*')
        else:
            if isinstance(dopr,Mv) and dopr.dot_flg:
                return Dop.Mul(self, dopr, op='*')
            else:
                result = Dop.Mul(self, dopr, op='*', diff=False)
                result.dot_flg = True
                return result

        raise ValueError('!!!!Incompatible dot flags in * __mul__!!!!')

    def __rmul__(self, dopl):  # dopl * self

        """
        Effect of overdot on Dop __rmul__

        Mv          *  Dop          = Dop
        Mv          *  Dop/Dot      = Dop
        Mv/Dot      *  Dop/Dot      = Mv
        """

        if self.dot_flg:  # self is dotted
            if not self.cmpflg:
                self.cmpflg = True
            if dopl.dot_flg:  # dopl is dotted
                return Dop.Mul(dopl, self, op='*')
            else:
                result = Dop.Mul(dopl, self, op='*', diff=False)
                result.dot_flg = True
                return result
        else:
            return Dop.Mul(dopl, self, op='*')

        raise ValueError('!!!!Incompatible dot flags in * __rmul__!!!!')

    def __truediv__(self, dopr):
        if isinstance(dopr, (Dop, Mv)):
            raise ValueError('In Dop.__truediv__ dopr must be a sympy scalar.')
        terms = []
        for term in self.terms:
            terms.append((term[0]/dopr,term[1]))
        return Dop(terms, ga= self.Ga)

    def __xor__(self, dopr):  # Dop (^) outer product

        if not self.dot_flg:
            return Dop.Mul(self, dopr, op='^')
        else:
            if isinstance(dopr,Mv) and dopr.dot_flg:
                return Dop.Mul(self, dopr, op='^')
            else:
                result = Dop.Mul(self, dopr, op='^', diff=False)
                result.dot_flg = True
                return result

        raise ValueError('!!!!Incompatible dot flags in ^ __xor__!!!!')

    def __rxor__(self, dopl):  # Dop (^) outer product

        if self.dot_flg:  # self is dotted
            if not self.cmpflg:
                self.cmpflg = True
            if dopl.dot_flg:  # dopl is dotted
                return Dop.Mul(dopl, self, op='^')
            else:
                result = Dop.Mul(dopl, self, op='^', diff=False)
                result.dot_flg = True
                return result
        else:
            return Dop.Mul(dopl, self, op='^')

        raise ValueError('!!!!Incompatible dot flags in ^ __rxor__!!!!')

    def __or__(self, dopr):  # | inner product

        if not self.dot_flg:
            return Dop.Mul(self, dopr, op='|')
        else:
            if isinstance(dopr,Mv) and dopr.dot_flg:
                return Dop.Mul(self, dopr, op='|')
            else:
                result = Dop.Mul(self, dopr, op='|', diff=False)
                result.dot_flg = True
                return result

        raise ValueError('!!!!Incompatible dot flags in | __or__!!!!')

    def __ror__(self, dopl):  # | inner product

        if self.dot_flg:  # self is dotted
            if not self.cmpflg:
                self.cmpflg = True
            if dopl.dot_flg:  # dopl is dotted
                return Dop.Mul(dopl, self, op='|')
            else:
                result = Dop.Mul(dopl, self, op='|', diff=False)
                result.dot_flg = True
                return result
        else:
            return Dop.Mul(dopl, self, op='|')

        raise ValueError('!!!!Incompatible dot flags in | __ror__!!!!')


    def __lt__(self, dopr):  # < left contraction
        return Dop.Mul(self, dopr, op='<')

    def __gt__(self, dopr):  # > right contraction
        return Dop.Mul(self, dopr, op='>')

    def __eq__(self, dop):
        if isinstance(dop, Dop):
            if self.Ga != dop.Ga:
                return False

            self = Sdop.consolidate_coefs(self)
            dop = Sdop.consolidate_coefs(dop)
            if len(self.terms) != len(dop.terms):
                return False
            if set(self.terms) != set(dop.terms):
                return False
            return True
        else:
            return False

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter

        return Printer().doprint(self)

    def __repr__(self):
        return str(self)

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        #if r'\begin{align*}' not in latex_str:
        if r'\begin{aligned}' not in latex_str:
            if self.title is None:
                latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
            else:
                latex_str = r'\begin{equation*} ' + self.title + ' = ' + latex_str + r' \end{equation*}'
        else:
            if self.title is not None:
                latex_str = latex_str.replace('&',' ' + self.title + ' =&',1)
        return latex_str

    def is_scalar(self):
        for x in self.terms:
            if isinstance(x[0], Mv) and not x[0].is_scalar():
                return False
        return True

    def is_zero(self):
        for x in self.terms:
            if x[0] != 0:
                return False
        return True

    def components(self):
        dop_lst = []
        for (sdop, base) in self.Dop_mv_expand():
            new_coefs = []
            new_pdiffs = []
            for (coef, pdiff) in sdop.terms:
                if pdiff in new_pdiffs:
                    index = new_pdiffs.index(pdiff)
                    new_coefs[index] += coef * base
                else:
                    new_pdiffs.append(pdiff)
                    new_coefs.append(coef * base)
            new_coefs = [Mv(x, ga=self.Ga) for x in new_coefs]
            terms = list(zip(new_coefs, new_pdiffs))
            dop_lst.append(Dop(terms, ga=self.Ga))
        return tuple(dop_lst)

    def Dop_mv_expand(self, modes=None):
        coefs = []
        bases = []
        self.consolidate_coefs()

        for (coef, pdiff) in self.terms:
            if isinstance(coef, Mv) and not coef.is_scalar():
                mv_terms = metric.linear_expand(coef.obj, mode=False)
                for (mv_coef, mv_base) in mv_terms:
                    if mv_base in bases:
                        index = bases.index(mv_base)
                        coefs[index] += Sdop([(mv_coef, pdiff)], ga=self.Ga)
                    else:
                        bases.append(mv_base)
                        coefs.append(Sdop([(mv_coef, pdiff)], ga=self.Ga))
            else:
                if isinstance(coef, Mv):
                    mv_coef = coef.obj
                else:
                    mv_coef = coef
                if S(1) in bases:
                    index = bases.index(S(1))
                    coefs[index] += Sdop([(mv_coef, pdiff)], ga=self.Ga)
                else:
                    bases.append(S(1))
                    coefs.append(Sdop([(mv_coef, pdiff)], ga=self.Ga))
        if modes is not None:
            for i in range(len(coefs)):
                coefs[i] = coefs[i].simplify(modes)
        terms = list(zip(coefs, bases))
        return sorted(terms, key=lambda x: self.Ga.blades_lst0.index(x[1]))


    def format_dop_coef(self,coef,base,first=False,latex_flg=False):

        coef_str = str(coef)
        base_str = str(base)
        add_flg  = isinstance(coef, Add)

        if base_str == '1':  # Case of 0-partial derivative
            if coef_str[0] == '-' or first:
                return coef_str
            else:
                return '+' + coef_str

        if coef_str == '1':
            if first:
                return base_str
            else:
                return '+' + base_str


        if coef_str == '-1':
            return '-' + base_str

        if latex_flg:
            if isinstance(coef, Add):
                if coef_str[0] == '-':
                    coef_str = r'-\left (' + str(-coef) + r'\right )'
                else:
                    coef_str = r'\left (' + coef_str + r'\right )'
        else:
            if isinstance(coef, Add):
                if coef_str[0] == '-':
                    coef_str = '-(' + str(-coef) + ')'
                else:
                    coef_str = '(' + coef_str + ')'

        if not first:
            if len(coef_str) > 0:
                if coef_str[0] != '-':
                    coef_str = '+' + coef_str
            else:
                coef_str = '+' + coef_str

        if latex_flg:
            return coef_str + base_str
        else:
            return coef_str + '*' + base_str










    def Dop_str(self):

        if len(self.terms) == 0:
            return ' 0 '

        self = self.consolidate_coefs()
        self.terms = sorted(self.terms,key = functools.cmp_to_key(Dop.compare_terms))
        s = ''
        first = True

        for (coef,base) in self.terms:

            coef = printer.format_coef(coef, first = first)
            first = False
            s += coef + str(base)

        return s


    def Dop_latex_str(self):

        if len(self.terms) == 0:
            return ' 0 '

        self = self.consolidate_coefs()
        self.terms = sorted(self.terms,key = functools.cmp_to_key(Dop.compare_terms))
        s = ''
        first = True

        if printer.GaLatexPrinter.dop_fmt == 2:
            if not printer.isinteractive():
                s += r'\renewcommand{\arraystretch}{2} '
            s = r'\begin{array}{l} '

        if len(self.terms) == 1:
            terms = self.terms
        else:
            terms = self.terms[:-1]

        for (coef,base) in terms:
            base_str = str(base)
            if base_str == '1':
                base_str = ''

            coef_str = printer.format_coef(coef, first, latex_flg = True)
            first = False
            if printer.GaLatexPrinter.dop_fmt == 1:
                s += coef_str + base_str
            else:
                s += coef_str + base_str + r' \\ '

        if len(self.terms) > 1:
            (coef,base) = self.terms[-1]
            coef_str = printer.format_coef(coef, first, latex_flg = True)
            s += coef_str + str(base)

        if printer.GaLatexPrinter.dop_fmt == 2:
            s += r' \end{array}'
            if not printer.isinteractive():
                s += r'\renewcommand{\arraystretch}{1}'

        return s

    def Fmt(self, fmt=1):
        """
        Set format for printing of multivector differential operators -

            fmt = 1 - One differential operator per line
            fmt = 2 - One partial derivative per line

        """

        if printer.GaLatexPrinter.latex_flg:
            printer.GaLatexPrinter.dop_fmt = fmt
            latex_str = printer.GaLatexPrinter.latex(self)
            printer.GaLatexPrinter.dop_fmt = 1
            return latex_str
        else:
            printer.GaPrinter.dop_fmt = fmt
            text_str = str(self)
            printer.GaPrinter.dop_fmt = 1
            return text_str

        return

    @staticmethod
    def basic(ga):
        r_basis = list(ga.r_basis)

        if not ga.is_ortho:
            r_basis = [x / ga.e_sq for x in r_basis]
        if ga.norm:
            r_basis = [x / e_norm for (x, e_norm) in zip(r_basis, ga.e_norm)]

        ga.lgrad = Dop(r_basis, ga.pdx, ga=ga)
        ga.rgrad = Dop(r_basis, ga.pdx, ga=ga, cmpflg=true)
        return ga.lgrad, ga.rgrad

################################# Alan Macdonald's additions #########################


def Nga(x, prec=5):
    if isinstance(x, Mv):
        Px = Mv(x, ga=x.Ga)
        Px.obj = Nsympy(x.obj, prec)
        return(Px)
    else:
        return(Nsympy(x, prec))


def printeigen(M):    # Print eigenvalues, multiplicities, eigenvectors of M.
    evects = M.eigenvects()
    for i in range(len(evects)):                   # i iterates over eigenvalues
        print(('Eigenvalue =', evects[i][0], '  Multiplicity =', evects[i][1], ' Eigenvectors:'))
        for j in range(len(evects[i][2])):         # j iterates over eigenvectors of a given eigenvalue
            result = '['
            for k in range(len(evects[i][2][j])):  # k iterates over coordinates of an eigenvector
                result += str(trigsimp(evects[i][2][j][k]).evalf(3))
                if k != len(evects[i][2][j]) - 1:
                    result += ', '
            result += '] '
            print(result)
            return


def printGS(M, norm=False):  # Print Gram-Schmidt output.
    from sympy import GramSchmidt
    global N
    N = GramSchmidt(M, norm)
    result = '[ '
    for i in range(len(N)):
        result += '['
        for j in range(len(N[0])):
            result += str(trigsimp(N[i][j]).evalf(3))
            if j != len(N[0]) - 1:
                result += ', '
        result += '] '
        if j != len(N[0]) - 1:
            result += ' '
    result += ']'
    print(result)
    return


def printrref(matrix, vars="xyzuvwrs"):   # Print rref of matrix with variables.
    rrefmatrix = matrix.rref()[0]
    rows, cols = rrefmatrix.shape
    if len(vars) < cols - 1:
        print('Not enough variables.')
        return
    for i in range(rows):
        result = ''
        for j in range(cols - 1):
            result += str(rrefmatrix[i, j]) + vars[j]
            if j != cols - 2:
                result += ' + '
        result += ' = ' + str(rrefmatrix[i, cols - 1])
        print(result)
        return

def com(A, B):
    return ga.Ga.com(A, B)

def correlation(u, v, dec=3):  # Compute the correlation coefficient of vectors u and v.
    rows, cols = u.shape
    uave = 0
    vave = 0
    for i in range(rows):
        uave += u[i]
        vave += v[i]
    uave = uave / rows
    vave = vave / rows
    ulocal = u[:, :]  # Matrix copy
    vlocal = v[:, :]
    for i in range(rows):
        ulocal[i] -= uave
        vlocal[i] -= vave
    return ulocal.dot(vlocal) / (ulocal.norm() * vlocal.norm()). evalf(dec)


def cross(v1, v2):
    if v1.is_vector() and v2.is_vector() and v1.Ga.name == v2.Ga.name and v1.Ga.n == 3:
        return -v1.Ga.I() * (v1 ^ v2)
    else:
        raise ValueError(str(v1) + ' and ' + str(v2) + ' not compatible for cross product.')


def dual(A):
    if isinstance(A, Mv):
        return A.dual()
    else:
        raise ValueError('A not a multivector in dual(A)')


def even(A):
    if not isinstance(A,Mv):
        raise ValueError('A = ' + str(A) + ' not a multivector in even(A).')
    return A.even()


def odd(A):
    if not isinstance(A,Mv):
        raise ValueError('A = ' + str(A) + ' not a multivector in even(A).')
    return A.odd()


def exp(A,hint='-'):
    if isinstance(A,Mv):
        return A.exp(hint)
    else:
        return sympy_exp(A)


def grade(A, r=0):
    if isinstance(A, Mv):
        return A.grade(r)
    else:
        raise ValueError('A not a multivector in grade(A,r)')


def inv(A):
    if not isinstance(A,Mv):
        raise ValueError('A = ' + str(A) + ' not a multivector in inv(A).')
    return A.inv()


def norm(A, hint='+'):
    if isinstance(A, Mv):
        return A.norm(hint=hint)
    else:
        raise ValueError('A not a multivector in norm(A)')


def norm2(A):
    if isinstance(A, Mv):
        return A.norm2()
    else:
        raise ValueError('A not a multivector in norm(A)')


def proj(B, A):  # Project on the blade B the multivector A
    if isinstance(A,Mv):
        return A.project_in_blade(B)
    else:
        raise ValueError('A not a multivector in proj(B,A)')


def rot(itheta, A, hint='-'):  # Rotate by the 2-blade itheta the multivector A
    if isinstance(A,Mv):
        return A.rotate_multivector(itheta, hint)
    else:
        raise ValueError('A not a multivector in rotate(A,itheta)')


def refl(B, A):  #  Project on the blade B the multivector A
    if isinstance(A,Mv):
        return A.reflect_in_blade(B)
    else:
        raise ValueError('A not a multivector in reflect(B,A)')


def rev(A):
    if isinstance(A, Mv):
        return A.rev()
    else:
        raise ValueError('A not a multivector in rev(A)')


def scalar(A):
    if not isinstance(A,Mv):
        raise ValueError('A = ' + str(A) + ' not a multivector in inv(A).')
    return A.scalar()

################################# MV class for backward compatibility ###################

class MV(Mv):

    @staticmethod
    def convert_metric(gstr):
        if gstr[0] is '[' and gstr[-1] is ']':
            gstr_lst = gstr[1:-1].split(',')
            g = []
            for x in gstr_lst:
                g.append(int(x))
            return g
        else:
            return gstr

    @staticmethod
    def setup(basis, metric=None, coords=None, rframe=False, debug=False, curv=(None,None)):

        if isinstance(metric,str):
            metric = MV.convert_metric(metric)
        if curv != (None,None):
            MV.GA = ga.Ga(basis, g=None, coords=coords, X=curv[0], debug=debug)
        else:
            MV.GA = ga.Ga(basis, g=metric, coords=coords, X=curv[0], debug=debug)
        MV.I = MV.GA.i
        MV.metric = MV.GA.g
        if coords is not None:
            (MV.grad,MV.rgrad) = MV.GA.grads()
            return list(MV.GA.mv()) + [MV.grad]
        else:
            return list(MV.GA.mv())


    def __init__(self, base, mvtype, fct=False, blade_rep=True):
        Mv.__init__(self, base, mvtype, f=fct, ga=MV.GA)

    def Fmt(self, fmt=1, title=None):
        print(Mv.Fmt(self, fmt=fmt, title=title))
        return

def ReciprocalFrame(basis, mode='norm'):

    GA = basis[0].Ga
    dim = len(basis)
    indexes = tuple(range(dim))
    index = [()]

    for i in indexes[-2:]:
        index.append(tuple(combinations(indexes, i + 1)))

    MFbasis = []

    for igrade in index[-2:]:
        grade = []
        for iblade in igrade:
            blade = Mv(1, 'scalar', ga=GA)
            for ibasis in iblade:
                blade ^= basis[ibasis]
            blade = blade.trigsimp()
            grade.append(blade)
        MFbasis.append(grade)
    E = MFbasis[-1][0]
    E_sq = trigsimp((E * E).scalar(),)

    duals = copy.copy(MFbasis[-2])

    duals.reverse()
    sgn = 1
    rbasis = []
    for dual in duals:
        recpv = (sgn * dual * E).trigsimp()
        rbasis.append(recpv)
        sgn = -sgn

    if mode != 'norm':
        rbasis.append(E_sq)
    else:
        for i in range(dim):
            rbasis[i] = rbasis[i] / E_sq

    return tuple(rbasis)

def compare_pdiffs(pdop1, pdop2):  # compare two Pdops

    if pdop1.order > pdop2.order:
        return 1
    if pdop1.order < pdop2.order:
        return -1

    keys1 = list(pdop1.pdiffs.keys())
    keys2 = list(pdop2.pdiffs.keys())
    lkeys1 = len(keys1)
    lkeys2 = len(keys2)

    if lkeys1 == lkeys2:
        s1 = ''.join([str(pdop1.Ga.coords.index(x)) for x in keys1])
        s2 = ''.join([str(pdop1.Ga.coords.index(x)) for x in keys2])
        if s1 < s2:
            return -1
        else:
            return 1
    else:
        if lkeys1 < lkeys2:
            return 1
        else:
            return -1

if __name__ == "__main__":
    pass
