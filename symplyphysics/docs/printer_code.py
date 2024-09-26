"""
Symplyphysics code printer
"""

from typing import Any

from sympy import S, Basic, Mul, Number, Pow, Rational, StrPrinter, sift, sympify, E
from sympy.core.mul import _keep_coeff
from sympy.printing.precedence import precedence
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
                return f"{self._print(S.One)} / {self.parenthesize(expr.base, prec, strict=False)}"

        e = self.parenthesize(expr.exp, prec, strict=False)
        if self.printmethod == "_sympyrepr" and expr.exp.is_Rational and expr.exp.q != 1:
            # the parenthesized exp should be '(Rational(a, b))' so strip parens,
            # but just check to be sure.
            if e.startswith("(Rational"):
                return f"{self.parenthesize(expr.base, prec, strict=False)}^{e[1:-1]}"
        return f"{self.parenthesize(expr.base, prec, strict=False)}^{e}"

    def _print_Mul(self, expr: Any) -> str:
        # pylint: disable=too-many-branches, too-many-statements
        prec = precedence(expr)

        # Check for unevaluated Mul. In this case we need to make sure the
        # identities are visible, multiple Rational factors are not combined
        # etc so we display in a straight-forward form that fully preserves all
        # args and their order.
        args = expr.args
        if args[0] is S.One or any(
                isinstance(a, Number) or
                a.is_Pow and all(ai.is_Integer for ai in a.args)
                for a in args[1:]):
            d, n = sift(args, lambda x:
                isinstance(x, Pow) and bool(x.exp.as_coeff_Mul()[0] < 0),
                binary=True)
            for i, di in enumerate(d):
                if di.exp.is_Number:
                    e = -di.exp
                else:
                    e = list(di.exp.args)
                    e[0] = -e[0]
                d[i] = Pow(di.base, e, evaluate=False) if e - 1 else di.base

            pre = []
            # don't parenthesize first factor if negative
            if n and not n[0].is_Add and n[0].could_extract_minus_sign():
                pre = [self._print(n.pop(0))]

            nfactors = pre + [self.parenthesize(a, prec, strict=False)
                for a in n]
            if not nfactors:
                nfactors = ["1"]

            # don't parenthesize first of denominator unless singleton
            if len(d) > 1 and d[0].could_extract_minus_sign():
                pre = [self._print(d.pop(0))]
            else:
                pre = []
            dfactors = pre + [self.parenthesize(a, prec, strict=False)
                for a in d]

            n = " * ".join(nfactors)
            d = " * ".join(dfactors)
            if len(dfactors) > 1:
                return f"{n} / ({d})"
            if dfactors:
                return f"{n} / {d}"
            return n

        c, e = expr.as_coeff_Mul()
        if c < 0:
            expr = _keep_coeff(-c, e)
            sign = "-"
        else:
            sign = ""

        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        pow_paren = []  # Will collect all pow with more than one base element and exp = -1

        # use make_args in case expr was something like -x -> x
        args = Mul.make_args(expr)

        # Gather args for numerator/denominator
        def apow(i: Pow) -> Basic:
            b, e = i.as_base_exp()
            eargs = list(Mul.make_args(e))
            if eargs[0] is S.NegativeOne:
                eargs = eargs[1:]
            else:
                eargs[0] = -eargs[0]
            eargs = [sympify(e) for e in eargs]
            mul_obj = Mul(*eargs)
            if isinstance(i, Pow):
                return i.func(b, mul_obj, evaluate=False)
            return i.func(mul_obj, evaluate=False)
        for item in args:
            if (item.is_commutative and
                    isinstance(item, Pow) and
                    bool(item.exp.as_coeff_Mul()[0] < 0)):
                if item.exp is not S.NegativeOne:
                    b.append(apow(item))
                else:
                    if (len(item.args[0].args) != 1 and
                            isinstance(item.base, (Mul, Pow))):
                        # To avoid situations like #14160
                        pow_paren.append(item)
                    b.append(item.base)
            elif item.is_Rational and item is not S.Infinity:
                if item.p != 1:
                    a.append(Rational(item.p))
                if item.q != 1:
                    b.append(Rational(item.q))
            else:
                a.append(item)

        a = a or [S.One]

        a_str = [self.parenthesize(x, prec, strict=False) for x in a]
        b_str = [self.parenthesize(x, prec, strict=False) for x in b]

        # To parenthesize Pow with exp = -1 and having more than one Symbol
        for item in pow_paren:
            if item.base in b:
                b_str[b.index(item.base)] = f"({b_str[b.index(item.base)]})"

        if not b:
            return sign + " * ".join(a_str)
        if len(b) == 1:
            return sign + " * ".join(a_str) + " / " + b_str[0]
        mul_str = " * ".join(b_str)
        return sign + " * ".join(a_str) + f" / ({mul_str})"

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
        value, base = expr.args
        str_value = self._print(value)
        
        if base == E:
            return f"log({str_value})"
        
        str_base = self._print(base)
        return f"log({str_value}, {str_base})"


def code_str(expr: Any, **settings: Any) -> str:
    return SymbolCodePrinter(settings).doprint(expr)
