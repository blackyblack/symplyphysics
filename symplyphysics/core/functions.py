from typing import Any
from sympy import Expr, log as sym_log, E, sqrt as sym_sqrt, Mul as SymMul


def log(expr: Expr, base: Expr = E, **kwargs: Any) -> Expr:
    kwargs.setdefault("evaluate", False)
    return sym_log(expr, base, **kwargs)


def sqrt(expr: Expr, **kwargs: Any) -> Expr:
    kwargs.setdefault("evaluate", False)
    return sym_sqrt(expr, **kwargs)


def mul(*exprs: Expr, **kwargs: Any) -> Expr:
    kwargs.setdefault("evaluate", False)
    return SymMul(*exprs, **kwargs)
