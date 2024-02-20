from __future__ import annotations
from typing import Any, Optional, Sequence, Self
from sympy import S, Symbol as SymSymbol, Expr, Equality
from sympy.physics.units import Dimension
from sympy.core.function import UndefinedFunction
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol, pretty_use_unicode
from .id_generator import next_id


class DimensionSymbol:
    _dimension: Dimension
    _display_name: str

    def __init__(self, display_name: str, dimension: Dimension = Dimension(S.One)) -> None:
        self._dimension = dimension
        self._display_name = display_name

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    @property
    def display_name(self) -> str:
        return self._display_name


class Symbol(DimensionSymbol, SymSymbol):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        display_name: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        **assumptions: Any) -> Self:
        name = next_name("SYM") if display_name is None else next_name(display_name)
        obj = SymSymbol.__new__(cls, name, **assumptions)
        return obj

    def __init__(self,
        display_name: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        **_assumptions: Any) -> None:
        display_name = self.name if display_name is None else display_name
        super().__init__(display_name, dimension)


class Function(DimensionSymbol, UndefinedFunction):

    name: str

    # NOTE: Self type cannot be used in a metaclass and 'mcs' is a metaclass here
    def __new__(mcs,
        display_name: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        **options: Any) -> Function:
        name = next_name("FUN") if display_name is None else next_name(display_name)
        obj = UndefinedFunction.__new__(mcs, name, **options)
        return obj

    def __init__(cls,
        display_name: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        **_options: Any) -> None:
        display_name = cls.name if display_name is None else display_name
        super().__init__(display_name, dimension)


# Symbol and Function have generated names, hence their display is not readable.
# Use custom implementation of the PrettyPrinter to convert real symbol names
# to user fiendly names.


class SymbolPrinter(PrettyPrinter):

    def __init__(self, **settings: Any) -> None:
        super().__init__(settings)

    def is_unicode(self) -> bool:
        return self._settings["use_unicode"]

    def _print_Symbol(self, e: Expr, bold_name: bool = False) -> prettyForm:
        symb_name = e.display_name if isinstance(e, Symbol) else getattr(e, "name")
        symb = pretty_symbol(symb_name, bold_name)
        return prettyForm(symb)

    # pylint: disable-next=too-many-arguments
    def _print_Function(self,
        e: Expr,
        sort: bool = False,
        func_name: Optional[str] = None,
        left: str = "(",
        right: str = ")") -> prettyForm:
        # optional argument func_name for supplying custom names
        # works only for applied functions
        func_name = e.func.display_name if isinstance(e.func, Function) else func_name
        return self._helper_print_function(e.func,
            e.args,
            sort=sort,
            func_name=func_name,
            left=left,
            right=right)

    # pylint: disable-next=invalid-name
    def _print_SumArray(self, e: Expr) -> prettyForm:
        return self._print_Function(e, func_name="SumArray")


def next_name(name: str) -> str:
    return name + str(next_id(name))


def print_expression(expr: Expr | Equality | Sequence[Expr | Equality]) -> str:
    pprinter = SymbolPrinter(use_unicode=False)
    # this is an ugly hack, but at least it works
    use_unicode = pprinter.is_unicode()
    uflag = pretty_use_unicode(use_unicode)
    try:
        return pprinter.doprint(expr)
    finally:
        pretty_use_unicode(uflag)


# Helper method for easier interaction with SumArray
def tuple_of_symbols(display_name: str,
    dimension: Dimension = Dimension(S.One),
    length: int = 1) -> tuple[Symbol, ...]:
    return tuple(Symbol(display_name + "_" + str(i), dimension) for i in range(length))
