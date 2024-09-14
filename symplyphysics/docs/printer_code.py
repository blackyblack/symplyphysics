"""
Symplyphysics code printer
"""

from typing import Any

from sympy import StrPrinter
from sympy.printing.precedence import precedence
from ..core.symbols.symbols import DimensionSymbolNew


class SymbolCodePrinter(StrPrinter):
    """
    A printer to convert Symplyphysics law expressions to symbols
    """

    def __init__(self, settings: Any = None) -> None:
        settings["order"] = "lex"
        StrPrinter.__init__(self, settings)

    # pylint: disable-next=invalid-name
    def _print_SymbolNew(self, expr: Any) -> str:
        return expr.display_name if isinstance(expr, DimensionSymbolNew) else getattr(expr, "name")

    def _print_Quantity(self, expr: Any) -> str:
        return self._print_SymbolNew(expr)

    # pylint: disable-next=invalid-name
    def _print_SymbolIndexedNew(self, expr: Any) -> str:
        return self._print_SymbolNew(expr)

    def _print_Pow(self, expr: Any, _rational: bool = False) -> str:
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
        if isinstance(expr, DimensionSymbolNew):
            return expr.display_name
        if isinstance(expr.func, DimensionSymbolNew):
            return expr.func.display_name
        args_str = self.stringify(expr.args, ", ")
        return expr.func.__name__ + f"({args_str})"

    # pylint: disable-next=invalid-name
    def _print_SumIndexed(self, expr: Any) -> str:
        # only one index of sum is supported
        # expr.args[0] contains indexed symbol with index applied
        # expr.args[0].args[0] contains just indexed symbol
        symbol, index = expr.args[0].args
        return f"Sum({self._print(symbol)}, {self._print(index)})"


def code_str(expr: Any, **settings: Any) -> str:
    return SymbolCodePrinter(settings).doprint(expr)
