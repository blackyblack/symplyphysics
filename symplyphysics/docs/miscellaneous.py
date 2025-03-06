from sympy import Expr, Mod, Integral, Product, Sum
from sympy.printing.precedence import precedence_traditional, PRECEDENCE


def needs_mul_brackets(expr: Expr, *, first: bool = False, last: bool = False) -> bool:
    """
    Checks whether ``expr`` needs to be wrapped in brackets when printed as part of `Mul`. This
    holds for `Add` and some container objects that would not need brackets when appearing last in
    `Mul`, e.g. `Integral`, `Product`, or `Sum`. The boolean flags ``first`` and ``last`` specify
    whether ``expr`` comes first or last in `Mul`, respectively.
    """

    if expr.is_Mul and not first and expr.could_extract_minus_sign():
        return True

    if precedence_traditional(expr) < PRECEDENCE["Mul"]:
        return True

    if expr.is_Relational or expr.is_Piecewise:
        return True

    if expr.has(Mod):
        return True

    if not last and expr.has(Integral, Product, Sum):
        return True

    return False


def needs_add_brackets(expr: Expr) -> bool:
    """
    Checks whether the expression needs to be wrapped in brackets when printed as part of `Add`.
    This is `False` for most things.
    """

    if expr.is_Relational:
        return True

    if expr.has(Mod):
        return True

    return False
