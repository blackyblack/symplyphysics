"""
Symplyphysics code printer
"""

from typing import Any

from sympy import StrPrinter
from sympy.printing.precedence import precedence
from ..core.symbols.symbols import DimensionSymbol


class SymbolCodePrinter(StrPrinter):
    """
    A printer to convert Symplyphysics law expressions to symbols
    """

    def __init__(self, settings: Any=None) -> None:
        StrPrinter.__init__(self, settings)

    def _print_Symbol(self, expr: Any) -> str:
        return expr.display_symbol if isinstance(expr, DimensionSymbol) else getattr(expr, "name")

    def _print_Quantity(self, expr: Any) -> str:
        return self._print_Symbol(expr)

    # pylint: disable-next=invalid-name
    def _print_SymbolIndexed(self, expr: Any) -> str:
        return self._print_Symbol(expr)

    def _print_Pow(self, expr: Any, _rational: bool=False) -> str:
        prec = precedence(expr)
        return f"{self.parenthesize(expr.base, prec)}^{self.parenthesize(expr.exp, prec)}"

    def _print_Relational(self, expr: Any) -> str:
        lhs_code = self._print(expr.lhs)
        rhs_code = self._print(expr.rhs)
        charmap = {
            "==": "=",
        }
        return f"{lhs_code} {charmap[expr.rel_op]} {rhs_code}"

    def _print_Function(self, expr: Any) -> str:
        if isinstance(expr, DimensionSymbol):
            return expr.display_symbol
        if isinstance(expr.func, DimensionSymbol):
            return expr.func.display_symbol
        return expr.func.__name__ + f"({self.stringify(expr.args, ", ")})"

    # pylint: disable-next=invalid-name
    def _print_SumIndexed(self, expr: Any) -> str:
        # only one index of sum is supported
        # expr.args[0] contains indexed symbol with index applied
        # expr.args[0].args[0] contains just indexed symbol
        return f"Sum({self._print(expr.args[0].args[0])}, i)"


def code_str(expr: Any, **settings: Any) -> str:
    return SymbolCodePrinter(settings).doprint(expr)
