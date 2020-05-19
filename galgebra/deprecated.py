import warnings

from . import ga
from .mv import Mv

# galgebra 0.5.0
warnings.warn(
    "The `galgebra.deprecated` module is deprecated",
    DeprecationWarning, stacklevel=2)

################################# MV class for backward compatibility ###################


class MV(Mv):
    """ A deprecated version of :class:`galgebra.mv.Mv`. """

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
    def setup(basis, metric=None, coords=None, rframe=False, debug=False, curv=(None, None)) -> list:
        """
        This function allows a single geometric algebra to be created.

        If the function is called more than once the old geometric algebra is
        overwritten by the new geometric algebra. The named input ``metric``
        is the same as the named input ``g`` in the current version of
        *galgebra*. Likewise, ``basis``, ``coords``, and ``debug`` are the same
        in the old and current versions of *galgebra* [17]_.
        Due to improvements in *sympy* the inputs ``rframe`` and ``curv[1]`` are
        no longer required. ``curv[0]`` is the vector function (list or tuple of
        scalar functions) of the coordinates required to define a vector manifold.
        For compatibility with the old version of *galgebra* if ``curv`` is used
        ``metric`` should be a orthonormal Euclidean metric of the same
        dimension as ``curv[0]``.

        It is strongly suggested that one use the new methods of defining a
        geometric algebra on a manifold.

        .. [17]
            If the metric is input as a list or list or lists the object is no
            longer quoted (input as a string). For example the old
            ``metric='[1,1,1]'`` becomes ``metric=[1,1,1]``.
        """

        if isinstance(metric, str):
            metric = MV.convert_metric(metric)
        if curv != (None, None):
            MV.GA = ga.Ga(basis, g=None, coords=coords, X=curv[0], debug=debug)
        else:
            MV.GA = ga.Ga(basis, g=metric, coords=coords, X=curv[0], debug=debug)
        MV.I = MV.GA.i
        MV.metric = MV.GA.g
        if coords is not None:
            MV.grad, MV.rgrad = MV.GA.grads()
            return list(MV.GA.mv()) + [MV.grad]
        else:
            return list(MV.GA.mv())

    def __init__(self, base, mvtype, fct=None, blade_rep=True):
        # galgebra 0.5.0
        warnings.warn(
            "The `galgebra.deprecated.MV` class is deprecated in favor of "
            "`galgebra.mv.Mv`.",
            DeprecationWarning, stacklevel=2)
        kwargs = {}
        if fct is not None:
            kwargs['f'] = fct  # only forward this argument if we received it
        Mv.__init__(self, base, mvtype, ga=MV.GA, **kwargs)

    def Fmt(self, fmt=1, title=None) -> None:
        """
        ``Fmt`` in ``MV`` has inputs identical to ``Fmt`` in ``Mv`` except that
        if ``A`` is a multivector then ``A.Fmt(2,'A')`` executes a print
        statement from ``MV`` and returns ``None``, while from ``Mv``,
        ``A.Fmt(2,'A')`` returns a string so that the function is compatible
        with use in *ipython notebook*.
        """
        print(Mv.Fmt(self, fmt=fmt, title=title))


def ReciprocalFrame(basis, mode='norm'):
    # galgebra 0.5.0
    warnings.warn(
        "The `galgebra.deprecated.ReciprocalFrame` function is deprecated in "
        "favor of the `ReciprocalFrame` method of `Ga` objects.",
        DeprecationWarning, stacklevel=2)
    GA = basis[0].Ga
    return GA.ReciprocalFrame(basis, mode=mode)
