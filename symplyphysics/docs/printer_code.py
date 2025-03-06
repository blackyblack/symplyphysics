"""
Symplyphysics code printer
"""

from typing import Any
from sympy import S, Expr, Mul, StrPrinter, E, Symbol as SymSymbol
from sympy.matrices.dense import DenseMatrix
from sympy.simplify import fraction
from sympy.printing.precedence import precedence
from ..core.symbols.symbols import DimensionSymbol, Function, IndexedSymbol
from .miscellaneous import needs_mul_brackets, needs_add_brackets


class SymbolCodePrinter(StrPrinter):  # type: ignore[misc]
    """
    A printer to convert Symplyphysics law expressions to symbols
    """

    def __init__(self, settings: Any = None) -> None:
        StrPrinter.__init__(self, settings)

    # pylint: disable-next=invalid-name
    def _print_Symbol(self, expr: Any) -> str:
        if isinstance(expr, DimensionSymbol):
            return expr.display_name

        return str(getattr(expr, "name"))

    # pylint: disable-next=invalid-name
    def _print_Quantity(self, expr: Any) -> str:
        return self._print_Symbol(expr)

    # pylint: disable-next=invalid-name
    def _print_IndexedSymbol(self, expr: Any) -> str:
        return self._print_Symbol(expr)

    # pylint: disable-next=invalid-name
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

    # pylint: disable-next=invalid-name
    def _print_Relational(self, expr: Any) -> str:
        lhs_code = self._print(expr.lhs)
        rhs_code = self._print(expr.rhs)
        charmap = {
            "==": "=",
        }
        return f"{lhs_code} {charmap[expr.rel_op]} {rhs_code}"

    # pylint: disable-next=invalid-name
    def _print_Function(self, expr: Any) -> str:
        if isinstance(expr, DimensionSymbol):
            if isinstance(expr, Function):
                name = expr.display_name
                args_str = self.stringify(expr.arguments, ", ")
            else:
                return expr.display_name
        else:
            func = expr.func
            if isinstance(func, DimensionSymbol):
                if isinstance(func, Function):
                    name = func.display_name
                else:
                    return func.display_name
            else:
                name = func.__name__
            args_str = self.stringify(expr.args, ", ")
        return f"{name}({args_str})"

    # pylint: disable-next=invalid-name
    def _print_IndexedSum(self, expr: Any) -> str:
        # only one index of sum is supported
        # expr.args[0] contains the argument of the Sum
        # expr.args[1] contains just indexed symbol
        applied, index = expr.args
        return f"Sum({self._print(applied)}, {self._print(index)})"

    # pylint: disable-next=invalid-name
    def _print_IndexedProduct(self, expr: Any) -> str:
        # only one index of product is supported
        # expr.args[0] contains the argument of the Product
        # expr.args[1] contains just indexed symbol
        applied, index = expr.args
        return f"Product({self._print(applied)}, {self._print(index)})"

    def _print_log(self, expr: Any) -> str:
        value, base = (expr.args[0], expr.args[1]) if len(expr.args) > 1 else (expr.args[0], E)
        str_value = self._print(value)

        if base == E:
            return f"log({str_value})"

        str_base = self._print(base)
        return f"log({str_value}, {str_base})"

    def _print_div(self, numer: Expr, denom: Expr) -> str:
        snumer = self._print_Mul(numer) if numer.is_Mul else str(self._print(numer))
        sdenom = self._print_Mul(denom) if denom.is_Mul else str(self._print(denom))

        snumer_str = f"({snumer})" if needs_mul_brackets(numer, first=True, last=False) else snumer
        mul_in_denom = False
        if denom.is_Mul:
            denom_args = [a for a in denom.args if a != S.One]
            mul_in_denom = len(denom_args) > 1
        sdenom_str = f"({sdenom})" if needs_mul_brackets(denom, first=False,
            last=True) or mul_in_denom else sdenom
        tex = f"{snumer_str} / {sdenom_str}"
        return tex

    # pylint: disable-next=invalid-name
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
                if needs_mul_brackets(term, first=i == 0, last=i == len(args) - 1):
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
                return tex + self._print_div(n_n, d)
            tex += convert_args(n_n)
            return tex
        tex += self._print_div(n_n, n_d)
        tex2 = ""
        if d_n == S.One:
            tex2 = convert_args(d_d)
        elif d_d == S.One:
            sdenom = convert_args(d_n)
            mul_in_denom = False
            if d_n.is_Mul:
                denom_args = [a for a in d_n.args if a != S.One]
                mul_in_denom = len(denom_args) > 1
            sdenom_str = (f"({sdenom})"
                if needs_mul_brackets(d_n, first=False, last=True) or mul_in_denom else sdenom)
            return f"{tex} / {sdenom_str}"
        else:
            tex2 = self._print_div(d_d, d_n)
        return f"{tex}{separator}{tex2}"

    # pylint: disable-next=invalid-name
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
            if needs_add_brackets(term):
                term_tex = f"({term_tex})"
            tex += term_tex

        return tex

    # pylint: disable-next=invalid-name
    def _print_DenseMatrix(self, expr: DenseMatrix) -> str:
        rows, cols = expr.shape

        if cols == 1 or rows == 1:
            parts = [self._print(elem) for elem in expr]
            result = "[" + ", ".join(parts) + "]"
            return result if cols == 1 else result + ".T"

        def print_row(row: int) -> str:
            parts = [self._print(expr[row, col]) for col in range(cols)]
            return "[" + ", ".join(parts) + "]"

        parts = [print_row(row) for row in range(rows)]
        return "[" + ", ".join(parts) + "]"

    # pylint: disable-next=invalid-name
    def _print_Average(self, expr: Any) -> str:
        return f"avg({self._print(expr.factor)})"

    # pylint: disable-next=invalid-name
    def _print_FiniteDifference(self, expr: Any) -> str:
        return f"Delta({self._print(expr.factor)})"

    # pylint: disable-next=invalid-name
    def _print_ExactDifferential(self, expr: Any) -> str:
        if expr.wrap_code:
            return f"d({self._print(expr.factor)})"

        return f"d{self._print(expr.factor)}"

    # pylint: disable-next=invalid-name
    def _print_InexactDifferential(self, expr: Any) -> str:
        return f"delta({self._print(expr.factor)})"


def code_str(expr: Any, **settings: Any) -> str:
    printer = SymbolCodePrinter(settings)

    if isinstance(expr, Function):
        arguments = [
            SymSymbol(arg.display_name) if isinstance(arg, Function) else arg
            for arg in expr.arguments
        ]
        expr = expr(*arguments)
    if isinstance(expr, IndexedSymbol):
        expr = expr[expr.index]

    return str(printer.doprint(expr))
