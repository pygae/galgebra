import copy
from itertools import combinations
from sympy import trigsimp, S
from . import ga
from .mv import Mv
from . import utils

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

        if utils.isstr(metric):
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


    def __init__(self, base, mvtype, fct=None, blade_rep=True):
        kwargs = {}
        if fct is not None:
            kwargs['f'] = fct  # only forward this argument if we received it
        Mv.__init__(self, base, mvtype, ga=MV.GA, **kwargs)

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
            blade = Mv(S(1), 'scalar', ga=GA)
            for ibasis in iblade:
                blade ^= basis[ibasis]
            blade = blade.trigsimp()
            grade.append(blade)
        MFbasis.append(grade)
    E = MFbasis[-1][0]
    E_sq = trigsimp((E * E).scalar(),)

    duals = copy.copy(MFbasis[-2])

    duals.reverse()
    sgn = S(1)
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