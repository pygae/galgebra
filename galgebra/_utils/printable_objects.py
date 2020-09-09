""" Small helper objects to aid in printing. """
import copy

from .printable import Printable


class WithSettings(Printable):
    """ Helper class to attach print settings to an object """
    def __init__(self, obj, settings: dict = {}):
        self._obj = obj
        self._settings = settings

    def __do_print(self, printer):
        # make a copy of the printer with the specified setting applied
        new_printer = copy.copy(printer)
        new_printer._settings = copy.copy(new_printer._settings)
        new_printer._settings.update(self._settings)
        return new_printer._print(self._obj)

    _latex = _pretty = _sympystr = __do_print


class FmtResult(Printable):
    """ Object returned from .Fmt methods, which can be printed as latex """
    def __new__(cls, obj: Printable, label: str) -> Printable:
        if label is None:
            return obj
        self = super().__new__(cls)
        self._obj = obj
        self._label = label
        return self

    def _latex(self, printer):
        return self._label + ' = ' + printer._print(self._obj)

    def _sympystr(self, printer):
        return self._label + ' = ' + printer._print(self._obj)
