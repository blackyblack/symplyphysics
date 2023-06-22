import random
import string
import sympy
from sympy.physics.units import Dimension
from sympy.core.function import UndefinedFunction
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol, pretty_use_unicode


class DimensionSymbol:
    _dimension: Dimension = None
    _display_name: str = None

    @staticmethod
    def random_name(base: str = None, digits: int = 4) -> str:
        return ("" if base is None else base) + "".join(random.choices(string.digits, k=digits))

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    @property
    def display_name(self) -> str:
        return self._display_name

    def __repr__(self) -> str:
        return self.display_name


class Symbol(DimensionSymbol, sympy.Symbol):

    def __new__(cls, display_name: str = None, dimension: Dimension = 1, **assumptions):
        name = DimensionSymbol.random_name("S",
            8) if display_name is None else DimensionSymbol.random_name(display_name)
        self = super().__new__(cls, name, **assumptions)
        self._dimension = dimension
        self._display_name = name if display_name is None else display_name
        return self


class Function(DimensionSymbol, UndefinedFunction):

    def __new__(cls, display_name: str = None, dimension: Dimension = 1, **options):
        name = DimensionSymbol.random_name("F",
            8) if display_name is None else DimensionSymbol.random_name(display_name)
        self = super().__new__(cls, name, **options)
        self._dimension = dimension
        self._display_name = name if display_name is None else display_name
        return self


# Symbol and Function have random name, hence their display is not readable.
# Use custom implementation of the PrettyPrinter to convert real symbol names
# to user fiendly names.


class SymbolPrinter(PrettyPrinter):

    def __init__(self, **settings):
        super().__init__(settings)

    def _print_Symbol(self, e, bold_name=False):
        symb_name = e.display_name if isinstance(e, Symbol) else e.name
        symb = pretty_symbol(symb_name, bold_name)
        return prettyForm(symb)

    def _print_Function(self, e, sort=False, func_name=None, left='(', right=')'):
        # optional argument func_name for supplying custom names
        # XXX works only for applied functions
        func_name = e.func.display_name if isinstance(e.func, Function) else func_name
        return self._helper_print_function(e.func,
            e.args,
            sort=sort,
            func_name=func_name,
            left=left,
            right=right)
    
    def _print_SumArray(self, expr):
        return self._print_Function(expr, func_name='SumArray')


def print_expression(expr) -> str:
    pp = SymbolPrinter(use_unicode=False)
    # XXX: this is an ugly hack, but at least it works
    use_unicode = pp._settings['use_unicode']
    uflag = pretty_use_unicode(use_unicode)
    try:
        return pp.doprint(expr)
    finally:
        pretty_use_unicode(uflag)


# Helper method for easier interaction with SumArray
def tuple_of_symbols(display_name: str = None, dimension: Dimension = 1, length: int = 1) -> bool:
    return tuple(Symbol(display_name + "_" + str(i), dimension) for i in range(length))
