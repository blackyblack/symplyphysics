from __future__ import annotations
from typing import Any, Optional, Sequence, Self
from sympy import S, Idx, Symbol as SymSymbol, Expr, Equality, IndexedBase
from sympy.physics.units import Dimension
from sympy.core.function import UndefinedFunction
from sympy.printing.printer import Printer
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol, pretty_use_unicode
from .id_generator import next_id


class DimensionSymbol:
    _dimension: Dimension
    _display_name: str
    _display_symbol: str
    _display_latex: str

    def __init__(self,
        display_name: str,
        dimension: Dimension = Dimension(S.One),
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None) -> None:
        self._dimension = dimension
        self._display_name = display_name
        self._display_symbol = self._display_name if display_symbol is None else display_symbol
        self._display_latex = self._display_symbol if display_latex is None else display_latex

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    @property
    def display_name(self) -> str:
        return self._display_name

    @property
    def display_symbol(self) -> str:
        return self._display_symbol

    @property
    def display_latex(self) -> str:
        return self._display_latex


class DimensionSymbolNew:
    _dimension: Dimension
    _display_name: str
    _display_latex: str

    def __init__(self,
        display_name: str,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None) -> None:
        self._dimension = dimension
        self._display_name = display_name
        self._display_latex = self._display_name if display_latex is None else display_latex

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    @property
    def display_name(self) -> str:
        return self._display_name

    @property
    def display_latex(self) -> str:
        return self._display_latex


class Symbol(DimensionSymbol, SymSymbol):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        display_name: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        **assumptions: Any) -> Self:
        name = next_name("SYM") if display_name is None else next_name(display_name)
        obj = SymSymbol.__new__(cls, name, **assumptions)
        return obj

    def __init__(self,
        display_name: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        **_assumptions: Any) -> None:
        display_name = self.name if display_name is None else display_name
        super().__init__(str(display_name),
            dimension,
            display_symbol=display_symbol,
            display_latex=display_latex)


class SymbolNew(DimensionSymbolNew, SymSymbol):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        display_symbol: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **assumptions: Any) -> Self:
        return SymSymbol.__new__(cls, next_name("SYM"), **assumptions)

    def __init__(self,
        display_symbol: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **_assumptions: Any) -> None:
        display_name = self.name if display_symbol is None else display_symbol
        super().__init__(str(display_name), dimension, display_latex=display_latex)


class SymbolIndexed(DimensionSymbol, IndexedBase):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        name_or_symbol: Optional[str | SymSymbol] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        **assumptions: Any) -> Self:
        # SymPy subs() and solve() creates dummy symbols. Allow create new indexed symbols
        # without renaming
        if isinstance(name_or_symbol, SymSymbol):
            return IndexedBase.__new__(cls, name_or_symbol, **assumptions)
        name = next_name("SYM") if name_or_symbol is None else next_name(name_or_symbol)
        obj = IndexedBase.__new__(cls, name, **assumptions)
        return obj

    def __init__(self,
        name_or_symbol: Optional[str | SymSymbol] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        **_assumptions: Any) -> None:
        display_name = self.name if name_or_symbol is None else str(name_or_symbol)
        super().__init__(display_name,
            dimension,
            display_symbol=display_symbol,
            display_latex=display_latex)

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return p.doprint(self.label if self.display_symbol is None else self.display_symbol)


class SymbolIndexedNew(DimensionSymbolNew, IndexedBase):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        name_or_symbol: Optional[str | SymSymbol] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **assumptions: Any) -> Self:
        # SymPy subs() and solve() creates dummy symbols. Allow create new indexed symbols
        # without renaming
        if isinstance(name_or_symbol, SymSymbol):
            return IndexedBase.__new__(cls, name_or_symbol, **assumptions)
        return IndexedBase.__new__(cls, next_name("SYM"), **assumptions)

    def __init__(self,
        name_or_symbol: Optional[str | SymSymbol] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **_assumptions: Any) -> None:
        display_name = str(self.name) if name_or_symbol is None else str(name_or_symbol)
        super().__init__(display_name, dimension, display_latex=display_latex)

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return p.doprint(self.label if self.display_name is None else self.display_name)


class Function(DimensionSymbol, UndefinedFunction):

    name: str

    # NOTE: Self type cannot be used in a metaclass and 'mcs' is a metaclass here
    def __new__(mcs,
        display_name: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        **options: Any) -> Function:
        name = next_name("FUN") if display_name is None else next_name(display_name)
        obj = UndefinedFunction.__new__(mcs, name, **options)
        return obj

    def __init__(cls,
        display_name: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        **_options: Any) -> None:
        display_name = cls.name if display_name is None else display_name
        super().__init__(str(display_name),
            dimension,
            display_symbol=display_symbol,
            display_latex=display_latex)


class FunctionNew(DimensionSymbolNew, UndefinedFunction):

    name: str

    # NOTE: Self type cannot be used in a metaclass and 'mcs' is a metaclass here
    def __new__(mcs,
        display_symbol: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **options: Any) -> Function:
        return UndefinedFunction.__new__(mcs, next_name("FUN"), **options)

    def __init__(cls,
        display_symbol: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **_options: Any) -> None:
        display_name = str(cls.name) if display_symbol is None else display_symbol
        super().__init__(display_name, dimension, display_latex=display_latex)


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

    # pylint: disable-next=invalid-name
    def _print_SymbolIndexed(self, e: Expr, bold_name: bool = False) -> prettyForm:
        symb_name = e.display_name if isinstance(e, SymbolIndexed) else getattr(e, "name")
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
    def _print_SumIndexed(self, e: Expr) -> prettyForm:
        return self._print_Function(e, func_name="SumIndexed")


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


def clone_symbol(source: SymbolNew,
    *,
    display_symbol: Optional[str] = None,
    display_latex: Optional[str] = None,
    **assumptions: Any) -> SymbolNew:
    assumptions = source.assumptions0 if assumptions is None or len(
        assumptions) == 0 else assumptions
    display_symbol = source.display_name if display_symbol is None else display_symbol
    display_latex = source.display_latex if display_latex is None else display_latex
    return SymbolNew(display_symbol, source.dimension, display_latex=display_latex, **assumptions)


# This is default index for indexed parameters, eg for using in SumIndexed
global_index = Idx("global_index")
