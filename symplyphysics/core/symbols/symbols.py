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

    def __init__(self,
        display_name: str,
        dimension: Dimension = Dimension(S.One)) -> None:
        self._dimension = dimension
        self._display_name = display_name

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    @property
    def display_name(self) -> str:
        return self._display_name


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

    def _sympystr(self, p: Printer) -> str:
        return p.doprint(self.display_name)


class Symbol(DimensionSymbol, SymSymbol):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        display_name: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        **assumptions: Any) -> Self:
        name = next_name("SYM") if display_name is None else next_name(display_name)
        return SymSymbol.__new__(cls, name, **assumptions)

    def __init__(self,
        display_name: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        **_assumptions: Any) -> None:
        display_name = str(self.name) if display_name is None else display_name
        super().__init__(display_name, dimension)


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
        display_name = str(self.name) if display_symbol is None else display_symbol
        super().__init__(display_name, dimension, display_latex=display_latex)


class SymbolIndexed(DimensionSymbol, IndexedBase):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        name_or_symbol: Optional[str | SymSymbol] = None,
        _dimension: Dimension = Dimension(S.One),
        **assumptions: Any) -> Self:
        # SymPy subs() and solve() creates dummy symbols. Allow create new indexed symbols
        # without renaming
        if isinstance(name_or_symbol, SymSymbol):
            return IndexedBase.__new__(cls, name_or_symbol, **assumptions)
        name = next_name("SYM") if name_or_symbol is None else next_name(name_or_symbol)
        return IndexedBase.__new__(cls, name, **assumptions)

    def __init__(self,
        name_or_symbol: Optional[str | SymSymbol] = None,
        dimension: Dimension = Dimension(S.One),
        **_assumptions: Any) -> None:
        display_name = str(self.name) if name_or_symbol is None else str(name_or_symbol)
        super().__init__(display_name, dimension)

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return p.doprint(self.display_name)


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


class Function(DimensionSymbol, UndefinedFunction):

    # NOTE: Self type cannot be used in a metaclass and 'mcs' is a metaclass here
    def __new__(mcs,
        display_name: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        **options: Any) -> Function:
        name = next_name("FUN") if display_name is None else next_name(display_name)
        return UndefinedFunction.__new__(mcs, name, **options)

    def __init__(cls,
        display_name: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        **_options: Any) -> None:
        display_name = str(cls.name) if display_name is None else display_name  # type: ignore[attr-defined]
        super().__init__(display_name, dimension)


class FunctionNew(DimensionSymbolNew, UndefinedFunction):

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
        display_name = str(cls.name) if display_symbol is None else display_symbol  # type: ignore[attr-defined]
        super().__init__(display_name, dimension, display_latex=display_latex)

    def __repr__(cls) -> str:
        return str(cls.display_name)


# Symbol and Function have generated names, hence their display is not readable.
# Use custom implementation of the PrettyPrinter to convert real symbol names
# to user fiendly names.


class SymbolPrinter(PrettyPrinter):

    def __init__(self, **settings: Any) -> None:
        super().__init__(settings)

    def is_unicode(self) -> bool:
        return self._settings["use_unicode"]

    def _print_Symbol(self, e: Expr, bold_name: bool = False) -> prettyForm:
        symb_name = e.display_name if isinstance(e, (Symbol, SymbolNew)) else getattr(e, "name")
        symb = pretty_symbol(symb_name, bold_name)
        return prettyForm(symb)

    # pylint: disable-next=invalid-name
    def _print_SymbolIndexed(self, e: Expr, bold_name: bool = False) -> prettyForm:
        return self._print_Symbol(e, bold_name)

    # pylint: disable-next=too-many-arguments
    def _print_Function(self,
        e: Expr,
        sort: bool = False,
        func_name: Optional[str] = None,
        left: str = "(",
        right: str = ")") -> prettyForm:
        # optional argument func_name for supplying custom names
        # works only for applied functions
        func_name = e.func.display_name if isinstance(e.func, (Function, FunctionNew)) else func_name
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
    display_symbol_new = source.display_name if display_symbol is None else display_symbol
    display_latex_new = display_latex
    if display_latex_new is None:
        display_latex_new = display_symbol
    if display_latex_new is None:
        display_latex_new = source.display_latex
    return SymbolNew(display_symbol_new, source.dimension, display_latex=display_latex_new, **assumptions)


# This is default index for indexed parameters, eg for using in SumIndexed
global_index = Idx("global_index")
