__all__ = ['Printable']

try:
    from sympy.printing.defaults import Printable
except ImportError:
    # sympy < 1.7
    from sympy import Basic

    class _Trick_can_print_latex(set):
        """
        A class that tricks the `_can_print_latex` function in sympy < 1.7 into
        returning True for our types, by subclassing a supported type.

        In sympy >= 1.7, it already returns True for our types.
        """
        def __init__(self, value):
            self.value = value

        def __iter__(self):
            return iter([])

        def _latex(self, printer):
            return printer._print(self.value)

    class Printable:
        """ Backport of `sympy.printing.defaults.Printable` """
        def __str__(self):
            from sympy.printing.str import sstr
            return sstr(self, order=None)

        __repr__ = __str__

        def _repr_latex_(self):
            f = None
            try:
                # This isn't perfect, in principle there could be multiple
                # active IPython's with different configurations.
                ip = get_ipython()
            except NameError:
                # Not in IPython
                pass
            else:
                # Reuse the printer that was customized by init_printing, if
                # present.
                f = ip.display_formatter.formatters['text/latex'].type_printers.get(Basic)

            if f is None:
                # no customizations present or ipython not running
                f = Basic._repr_latex_

            if f is None:
                # latex printing disabled
                return None

            return f(_Trick_can_print_latex(self))
