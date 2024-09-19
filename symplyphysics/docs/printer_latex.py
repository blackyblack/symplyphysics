"""
Symplyphysics latex printer
"""

from typing import Any
from sympy.printing.latex import LatexPrinter, accepted_latex_functions
from sympy.core.function import AppliedUndef
from ..core.symbols.symbols import DimensionSymbolNew


class SymbolLatexPrinter(LatexPrinter):
    """
    A printer to convert Symplyphysics law expressions to latex
    """
    language = "Symplyphysics"

    def __init__(self, settings: Any = None) -> None:
        LatexPrinter.__init__(self, settings)

    # pylint: disable-next=invalid-name
    def _print_SymbolNew(self, expr: Any, style: str = "plain") -> str:
        display_name = expr.display_latex if isinstance(expr, DimensionSymbolNew) else getattr(
            expr, "name")
        name: str = self._settings["symbol_names"].get(display_name)
        if name is not None:
            return name
        return self._deal_with_super_sub(display_name, style=style)

    # pylint: disable-next=invalid-name
    def _print_Quantity(self, expr: Any) -> str:
        return self._print_SymbolNew(expr)

    # pylint: disable-next=invalid-name
    def _print_SymbolIndexedNew(self, expr: Any) -> str:
        return self._print_SymbolNew(expr)

    #pylint: disable-next=too-many-branches
    def _print_Function(self, expr: Any, exp: Any = None) -> str:
        func = expr.func.__name__
        if isinstance(expr, DimensionSymbolNew):
            func = expr.display_latex
        if isinstance(expr.func, DimensionSymbolNew):
            func = expr.func.display_latex

        if hasattr(self, "_print_" + func) and not isinstance(expr, AppliedUndef):
            return getattr(self, "_print_" + func)(expr, exp)

        args = [str(self._print(arg)) for arg in expr.args]
        # How inverse trig functions should be displayed, formats are:
        # abbreviated: asin, full: arcsin, power: sin^-1
        inv_trig_style = self._settings["inv_trig_style"]
        # If we are dealing with a power-style inverse trig function
        inv_trig_power_case = False
        # If it is applicable to fold the argument brackets
        can_fold_brackets = self._settings["fold_func_brackets"] and \
            len(args) == 1 and \
            not self._needs_function_brackets(expr.args[0])

        inv_trig_table = [
            "asin",
            "acos",
            "atan",
            "acsc",
            "asec",
            "acot",
            "asinh",
            "acosh",
            "atanh",
            "acsch",
            "asech",
            "acoth",
        ]

        # If the function is an inverse trig function, handle the style
        if func in inv_trig_table:
            if inv_trig_style == "abbreviated":
                pass
            elif inv_trig_style == "full":
                func = ("ar" if func[-1] == "h" else "arc") + func[1:]
            elif inv_trig_style == "power":
                func = func[1:]
                inv_trig_power_case = True

                # Can never fold brackets if we're raised to a power
                if exp is not None:
                    can_fold_brackets = False

        if inv_trig_power_case:
            if func in accepted_latex_functions:
                name = r"\%s^{-1}" % func
            else:
                name = r"\operatorname{%s}^{-1}" % func
        elif exp is not None:
            func_tex = self._hprint_Function(func)
            func_tex = self.parenthesize_super(func_tex)
            name = r"%s^{%s}" % (func_tex, exp)
        else:
            name = self._hprint_Function(func)

        if can_fold_brackets:
            if func in accepted_latex_functions:
                # Wrap argument safely to avoid parse-time conflicts
                # with the function name itself
                name += r" {%s}"
            else:
                name += r"%s"
        else:
            name += r"{\left(%s \right)}"

        if inv_trig_power_case and exp is not None:
            name += r"^{%s}" % exp

        return name % ",".join(args)

    # TODO: use e^ for shorter expressions
    def _print_ExpBase(self, expr: Any, exp: Any = None) -> str:
        args = [str(self._print(arg)) for arg in expr.args]
        can_fold_brackets = self._settings["fold_func_brackets"] and \
            len(args) == 1 and \
            not self._needs_function_brackets(expr.args[0])
        args_str = self._print(expr.args[0])
        name = f"{args_str}" if can_fold_brackets else f"\\left({args_str} \\right)"
        tex = r"\exp{%s}" % name
        return self._do_exponent(tex, exp)

    # pylint: disable-next=invalid-name
    def _print_SumIndexed(self, expr: Any) -> str:
        # only one index of sum is supported
        # expr.args[0] contains indexed symbol with index applied
        # expr.args[0].args[0] contains just indexed symbol
        symbol, index = expr.args[0].args
        return rf"\sum_{self._print(index)} {self._print(symbol)}"


def latex_str(expr: Any, **settings: Any) -> str:
    return SymbolLatexPrinter(settings).doprint(expr)
