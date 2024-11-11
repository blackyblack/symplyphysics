"""
Symplyphysics latex printer
"""

import re
from typing import Any
from sympy import E, S, Expr, Mod, Mul
from sympy.printing.latex import LatexPrinter, accepted_latex_functions
from sympy.core.function import AppliedUndef
from sympy.simplify import fraction
from ..core.symbols.symbols import DimensionSymbolNew, FunctionNew

_between_two_numbers_p = (
    re.compile(r"[0-9][} ]*$"),  # search
    re.compile(r"(\d|\\frac{\d+}{\d+})"),  # match
)


def _discard_minus_sign(expr: Expr) -> tuple[Expr, bool]:
    if not isinstance(expr, Mul):
        if expr.could_extract_minus_sign():
            return (-expr, True)
        return (expr, False)
    args = []
    sign = False
    for a in expr.args:
        a, arg_sign = _discard_minus_sign(a)
        args.append(a)
        if arg_sign:
            sign = not sign
    return (Mul(*args, evaluate=False), sign)


class SymbolLatexPrinter(LatexPrinter):
    """
    A printer to convert Symplyphysics law expressions to latex
    """
    language = "Symplyphysics"

    def __init__(self, settings: Any = None) -> None:
        settings["order"] = "none"
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
        index = expr.index
        str_expr = expr.display_latex
        str_index = self._print(index)
        return f"{{{str_expr}}}_{{{str_index}}}"

    def _print_Function(self, expr: Any, exp: Any = None) -> str:
        # pylint: disable=too-many-branches
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
        tex = f"\\exp{{{name}}}"
        return self._do_exponent(tex, exp)

    # pylint: disable-next=invalid-name
    def _print_SumIndexed(self, expr: Any) -> str:
        # only one index of sum is supported
        # expr.args[0] contains indexed symbol with index applied
        # expr.args[0].args[0] contains just indexed symbol
        symbol, index = expr.args[0].args
        return f"\\sum_{self._print(index)} {self._print(symbol)}"

    def _print_log(self, expr: Any, _exp: Any = None) -> str:
        value, base = (expr.args[0], expr.args[1]) if len(expr.args) > 1 else (expr.args[0], E)
        str_value = self._print(value)
        str_base = self._print(base)
        head = "\\log" if (base is None or base == E) else f"\\log_{{{str_base}}}"
        return f"{head} \\left( {str_value} \\right)"

    def _print_div(self, numer: Expr, denom: Expr) -> str:
        snumer = self._print_Mul(numer) if numer.is_Mul else str(self._print(numer))
        sdenom = self._print_Mul(denom) if denom.is_Mul else str(self._print(denom))
        tex = f"\\frac{{{snumer}}}{{{sdenom}}}"
        return tex

    def _extract_minus_sign(self, expr: Expr) -> tuple[Expr, bool]:
        args = expr.args
        sign = False
        if expr.is_Number and expr.is_extended_negative:
            return (-expr, True)
        if not expr.is_Mul:
            return (expr, False)
        terms = []
        for a in args:
            term, term_sign = self._extract_minus_sign(a)
            if term_sign:
                sign = not sign
            terms.append(term)
        return (Mul(*terms, evaluate=False), sign)


    def _print_Mul(self, expr: Mul) -> str:
        separator: str = self._settings["mul_symbol_latex"]
        numbersep: str = self._settings["mul_symbol_latex_numbers"]

        def convert_args(expr: Expr) -> str:
            if not expr.is_Mul:
                return str(self._print(expr))

            args = expr.args
            _tex = last_term_tex = ""

            if len(args) == 1 and args[0] == S.One:
                return "1"
            # Filter all 1 multiplications
            args = [a for a in args if a != S.One]

            for i, term in enumerate(args):
                term_tex = self._print(term)
                if self._needs_mul_brackets(term, first=i == 0, last=i == len(args) - 1):
                    term_tex = f"\\left({term_tex}\\right)"

                if  _between_two_numbers_p[0].search(last_term_tex) and \
                    _between_two_numbers_p[1].match(term_tex):
                    # between two numbers
                    _tex += numbersep
                elif _tex:
                    _tex += separator

                _tex += term_tex
                last_term_tex = term_tex
            return _tex

        expr, sign = self._extract_minus_sign(expr)
        tex = "- " if sign else ""

        n, d = fraction(expr, exact=True)

        n_n, n_d = fraction(n, exact=True)
        d_n, d_d = fraction(d, exact=True)

        # double fraction
        if S.One not in (n_d, d):
            tex += convert_args(expr)
            return tex

        # no denominator
        if n_d == S.One:
            if d != S.One:
                tex += self._print_div(n_n, d)
                return tex
            tex += convert_args(n_n)
            return tex
        tex = self._print_div(n_n, n_d)
        tex2 = ""
        if d_n == S.One:
            tex2 = convert_args(d_d)
        else:
            tex2 = self._print_div(d_d, d_n)
        return tex + separator + tex2

    def _needs_add_brackets(self, expr: Expr) -> bool:
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed as part of an Add, False otherwise.  This is False for most
        things.
        """
        if expr.is_Relational:
            return True
        if any(expr.has(x) for x in (Mod,)):
            return True
        return False

    def _print_Add(self, expr: Expr, _order: bool = False) -> str:
        tex = ""
        for i, term in enumerate(expr.args):
            if i == 0:
                pass
            elif term.could_extract_minus_sign():
                tex += " - "
                term = -term
            else:
                tex += " + "
            term_tex = self._print(term)
            if self._needs_add_brackets(term):
                term_tex = f"\\left({term_tex}\\right)"
            tex += term_tex

        return tex


def latex_str(expr: Any, **settings: Any) -> str:
    printer = SymbolLatexPrinter(settings)

    if isinstance(expr, FunctionNew):
        expr = expr(*expr.arguments)

    return printer.doprint(expr)
