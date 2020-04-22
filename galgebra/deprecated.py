from . import ga
from .mv import Mv

################################# MV class for backward compatibility ###################

class MV(Mv):

    @staticmethod
    def convert_metric(gstr):
        if gstr[0] == '[' and gstr[-1] == ']':
            gstr_lst = gstr[1:-1].split(',')
            g = []
            for x in gstr_lst:
                g.append(int(x))
            return g
        else:
            return gstr

    @staticmethod
    def setup(basis, metric=None, coords=None, rframe=False, debug=False, curv=(None,None)):

        if isinstance(metric, str):
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


def ReciprocalFrame(basis, mode='norm'):
    GA = basis[0].Ga
    return GA.ReciprocalFrame(basis, mode=mode)
