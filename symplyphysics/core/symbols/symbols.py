import sympy
from sympy.physics.units import Dimension
from sympy.core.function import UndefinedFunction
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol, pretty_use_unicode
from .id_generator import next_id


class DimensionSymbol:
    _dimension: Dimension = None
    _display_name: str = None

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
        name = next_name("S") if display_name is None else next_name(display_name)
        self = super().__new__(cls, name, **assumptions)
        self._dimension = dimension
        self._display_name = name if display_name is None else display_name
        return self


class Function(DimensionSymbol, UndefinedFunction):

    def __new__(cls, display_name: str = None, dimension: Dimension = 1, **options):
        name = next_name("F") if display_name is None else next_name(display_name)
        self = super().__new__(cls, name, **options)
        self._dimension = dimension
        self._display_name = name if display_name is None else display_name
        return self


# Symbol and Function have generated names, hence their display is not readable.
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


def next_name(name: str = None) -> str:
    return name + str(next_id(name))


def print_expression(expr) -> str:
    pp = SymbolPrinter(use_unicode=False)
    # XXX: this is an ugly hack, but at least it works
    use_unicode = pp._settings['use_unicode']
    uflag = pretty_use_unicode(use_unicode)
    try:
        return pp.doprint(expr)
    finally:
        pretty_use_unicode(uflag)
