from __future__ import annotations
from typing import Any, Optional, Sequence
from sympy import S, Idx, MatAdd, MatMul, MatrixBase, Symbol as SymSymbol, Expr, Equality, IndexedBase, Matrix as SymMatrix
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
    _display_latex: str

    def __init__(self,
        display_name: str,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None) -> None:
        self._dimension = dimension
        self._display_name = display_name
        self._display_latex = display_latex or self._display_name

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
        return str(p.doprint(self.display_name))


class Symbol(DimensionSymbol, SymSymbol):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        display_symbol: Optional[str] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **assumptions: Any) -> Symbol:
        obj = SymSymbol.__new__(cls, next_name("SYM"), **assumptions)
        return obj

    def __init__(self,
        display_symbol: Optional[str] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **_assumptions: Any) -> None:
        display_name = display_symbol or str(self.name)
        super().__init__(display_name, dimension, display_latex=display_latex)

    # HACK: fix for pylint false positive on unary minus or plus
    def __neg__(self) -> Symbol:  # pylint: disable=useless-parent-delegation
        return super().__neg__()


# This is default index for indexed parameters, e.g. for using in IndexedSum
global_index = Idx("i")


class IndexedSymbol(DimensionSymbol, IndexedBase):  # pylint: disable=too-many-ancestors
    index: Idx

    def __new__(cls,
        name_or_symbol: Optional[str | SymSymbol] = None,
        index: Optional[Idx] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **assumptions: Any) -> IndexedSymbol:
        # SymPy subs() and solve() creates dummy symbols. Allow create new indexed symbols
        # without renaming
        if isinstance(name_or_symbol, SymSymbol):
            obj = IndexedBase.__new__(cls, name_or_symbol, **assumptions)
        else:
            obj = IndexedBase.__new__(cls, next_name("SYM"), **assumptions)
        return obj

    def __init__(self,
        name_or_symbol: Optional[str | SymSymbol] = None,
        index: Optional[Idx] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **_assumptions: Any) -> None:
        display_name = str(self.name) if name_or_symbol is None else str(name_or_symbol)
        self.index = index or global_index
        super().__init__(display_name, dimension, display_latex=display_latex)

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass


class Function(DimensionSymbol, UndefinedFunction):
    arguments: Optional[Sequence[Expr]]

    # NOTE: Self type cannot be used in a metaclass and 'mcs' is a metaclass here
    # NOTE: constructor returns not an object, but a class. Object is constructed
    #       when arguments of a function are applied.
    def __new__(mcs,
        display_symbol: Optional[str] = None,
        arguments: Optional[Sequence[Expr]] = None,
        _dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **options: Any) -> Function:
        obj = UndefinedFunction.__new__(mcs, next_name("FUN"), **options)
        return obj

    def __init__(cls,
        display_symbol: Optional[str] = None,
        arguments: Optional[Sequence[Expr]] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
        **options: Any) -> None:
        display_name = display_symbol or str(cls.name)
        cls.arguments = arguments
        DimensionSymbol.__init__(cls, display_name, dimension, display_latex=display_latex)

        if arguments is not None:
            options["nargs"] = len(arguments)
        UndefinedFunction.__init__(cls, **options)

    def __repr__(cls) -> str:  # pylint: disable=invalid-repr-returned
        return cls.display_name


class Matrix(SymMatrix):  # pylint: disable=too-many-ancestors

    def __mul__(self: MatrixBase, other: MatrixBase) -> Expr:
        return MatMul(self, other)

    def __add__(self: MatrixBase, other: MatrixBase) -> Expr:
        return MatAdd(self, other)


# Symbol and Function have generated names, hence their display is not readable.
# Use custom implementation of the PrettyPrinter to convert real symbol names
# to user friendly names.


class SymbolPrinter(PrettyPrinter):

    def __init__(self, **settings: Any) -> None:
        super().__init__(settings)

    def is_unicode(self) -> bool:
        return bool(self._settings["use_unicode"])

    def _print_Symbol(self, e: Expr, bold_name: bool = False) -> prettyForm:
        symb_name = e.display_name if isinstance(e, Symbol) else getattr(e, "name")
        symb = pretty_symbol(symb_name, bold_name)
        return prettyForm(symb)

    # pylint: disable-next=invalid-name
    def _print_SymbolIndexed(self, e: Expr, bold_name: bool = False) -> prettyForm:
        return self._print_Symbol(e, bold_name)

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
    def _print_IndexedSum(self, e: Expr) -> prettyForm:
        return self._print_Function(e, func_name="IndexedSum")


def next_name(name: str) -> str:
    return name + str(next_id(name))


def print_expression(expr: Expr | Equality | Sequence[Expr | Equality]) -> str:
    pprinter = SymbolPrinter(use_unicode=False)
    # this is an ugly hack, but at least it works
    use_unicode = pprinter.is_unicode()
    uflag = pretty_use_unicode(use_unicode)
    try:
        return str(pprinter.doprint(expr))
    finally:
        pretty_use_unicode(uflag)


def process_subscript_and_names(
    code_name: str,
    latex_name: str,
    subscript: Optional[str] = None,
) -> tuple[str, str]:
    if not subscript:
        return code_name, latex_name

    return f"{code_name}_{subscript}", f"{latex_name}_{{{subscript}}}"


def clone_as_symbol(source: Symbol | IndexedSymbol,
    *,
    display_symbol: Optional[str] = None,
    display_latex: Optional[str] = None,
    subscript: Optional[str] = None,
    **assumptions: Any) -> Symbol:
    assumptions = assumptions or source.assumptions0
    display_symbol = display_symbol or source.display_name
    display_latex = display_latex or source.display_latex

    display_symbol, display_latex = process_subscript_and_names(display_symbol, display_latex,
        subscript)

    return Symbol(
        display_symbol,
        source.dimension,
        display_latex=display_latex,
        **assumptions,
    )


def clone_as_function(
    source: Symbol | IndexedSymbol,
    arguments: Optional[Sequence[Expr]] = None,
    *,
    display_symbol: Optional[str] = None,
    display_latex: Optional[str] = None,
    subscript: Optional[str] = None,
    **assumptions: Any,
) -> Function:
    display_symbol = display_symbol or source.display_name
    display_latex = display_latex or source.display_latex

    display_symbol, display_latex = process_subscript_and_names(display_symbol, display_latex,
        subscript)

    return Function(
        display_symbol,
        arguments,
        source.dimension,
        display_latex=display_latex,
        **assumptions,
    )


def clone_as_indexed(
    source: Symbol | IndexedSymbol,
    index: Optional[Idx] = None,
    *,
    display_symbol: Optional[str] = None,
    display_latex: Optional[str] = None,
    **assumptions: Any,
) -> IndexedSymbol:
    assumptions = assumptions or source.assumptions0
    display_symbol = display_symbol or source.display_name
    display_latex = display_latex or source.display_latex
    return IndexedSymbol(
        display_symbol,
        index,
        source.dimension,
        display_latex=display_latex,
        **assumptions,
    )
