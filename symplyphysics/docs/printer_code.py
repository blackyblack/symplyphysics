"""
Symplyphysics code printer
"""

from typing import Any
from sympy import S, Expr, Mod, Mul, StrPrinter, E
from sympy.simplify import fraction
from sympy.concrete.products import Product
from sympy.concrete.summations import Sum
from sympy.integrals.integrals import Integral
from sympy.printing.precedence import precedence, precedence_traditional, PRECEDENCE
from ..core.symbols.symbols import DimensionSymbolNew


class SymbolCodePrinter(StrPrinter):
    """
    A printer to convert Symplyphysics law expressions to symbols
    """

    def __init__(self, settings: Any = None) -> None:
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

        if expr.exp is S.Half:
            return f"sqrt({self._print(expr.base)})"

        if expr.is_commutative:
            if -expr.exp is S.Half:
                # Note: Don't test "expr.exp == -S.Half" here, because that will
                # match -0.5, which we don't want.
                n, d = (self._print(arg) for arg in (S.One, expr.base))
                return f"{n} / sqrt({d})"
            if expr.exp is -S.One:
                # Similarly to the S.Half case, don't test with "==" here.
                return f"{self._print(S.One)} / {self.parenthesize(expr.base, prec, strict=True)}"

        e = self.parenthesize(expr.exp, prec, strict=False)
        if self.printmethod == "_sympyrepr" and expr.exp.is_Rational and expr.exp.q != 1:
            # the parenthesized exp should be '(Rational(a, b))' so strip parens,
            # but just check to be sure.
            if e.startswith("(Rational"):
                return f"{self.parenthesize(expr.base, prec, strict=False)}^{e[1:-1]}"
        return f"{self.parenthesize(expr.base, prec, strict=False)}^{e}"

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

    def _print_log(self, expr: Any) -> str:
        value, base = (expr.args[0], expr.args[1]) if len(expr.args) > 1 else (expr.args[0], E)
        str_value = self._print(value)

        if base == E:
            return f"log({str_value})"

        str_base = self._print(base)
        return f"log({str_value}, {str_base})"

    def _needs_mul_brackets(self, expr: Expr, first: bool = False, last: bool = False) -> bool:
        # pylint: disable=too-many-return-statements
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed as part of a Mul, False otherwise. This is True for Add,
        but also for some container objects that would not need brackets
        when appearing last in a Mul, e.g. an Integral. ``last=True``
        specifies that this expr is the last to appear in a Mul.
        ``first=True`` specifies that this expr is the first to appear in
        a Mul.
        """
        if expr.is_Mul:
            if not first and expr.could_extract_minus_sign():
                return True
        elif precedence_traditional(expr) < PRECEDENCE["Mul"]:
            return True
        elif expr.is_Relational:
            return True
        if expr.is_Piecewise:
            return True
        if any(expr.has(x) for x in (Mod,)):
            return True
        if (not last and any(expr.has(x) for x in (Integral, Product, Sum))):
            return True

        return False

    def _print_div(self, numer: Expr, denom: Expr) -> str:
        snumer = self._print_Mul(numer) if numer.is_Mul else str(self._print(numer))
        sdenom = self._print_Mul(denom) if denom.is_Mul else str(self._print(denom))

        snumer_str = f"({snumer})" if self._needs_mul_brackets(numer, True, False) else snumer
        mul_in_denom = False
        if denom.is_Mul:
            denom_args = [a for a in denom.args if a != S.One]
            mul_in_denom = len(denom_args) > 1
        sdenom_str = f"({sdenom})" if self._needs_mul_brackets(denom, False,
            True) or mul_in_denom else sdenom
        tex = f"{snumer_str} / {sdenom_str}"
        return tex

    def _print_Mul(self, expr: Mul) -> str:
        separator = " * "

        def convert_args(expr: Expr) -> str:
            if not expr.is_Mul:
                return str(self._print(expr))

            args = expr.args
            _tex = ""

            if len(args) == 1 and args[0] == S.One:
                return "1"
            # Filter all 1 multiplications
            args = [a for a in args if a != S.One]

            for i, term in enumerate(args):
                term_tex = self._print(term)
                if self._needs_mul_brackets(term, first=i == 0, last=i == len(args) - 1):
                    term_tex = f"({term_tex})"

                if _tex:
                    _tex += separator

                _tex += term_tex
            return _tex

        tex = ""
        if expr.could_extract_minus_sign():
            expr = -expr
            tex = "-"

        n, d = fraction(expr, exact=True)

        n_n, n_d = fraction(n, exact=True)
        d_n, d_d = fraction(d, exact=True)

        # no denominator
        if n_d == S.One:
            if d != S.One:
                return self._print_div(n_n, d)
            tex += convert_args(n_n)
            return tex
        tex = self._print_div(n_n, n_d)
        tex2 = ""
        if d_n == S.One:
            tex2 = convert_args(d_d)
        elif d_d == S.One:
            sdenom = convert_args(d_n)
            mul_in_denom = False
            if d_n.is_Mul:
                denom_args = [a for a in d_n.args if a != S.One]
                mul_in_denom = len(denom_args) > 1
            sdenom_str = f"({sdenom})" if self._needs_mul_brackets(d_n, False,
                True) or mul_in_denom else sdenom
            return f"{tex} / {sdenom_str}"
        else:
            tex2 = self._print_div(d_d, d_n)
        return f"{tex}{separator}{tex2}"

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
                term_tex = f"({term_tex})"
            tex += term_tex

        return tex


def code_str(expr: Any, **settings: Any) -> str:
    return SymbolCodePrinter(settings).doprint(expr)
