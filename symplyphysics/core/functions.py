from sympy import Expr, log as sym_log, E, sqrt as sym_sqrt, Mul


def log(expr: Expr, base: Expr = E, **kwargs) -> Expr:
    kwargs.setdefault("evaluate", False)
    return sym_log(expr, base, **kwargs)


def sqrt(expr: Expr, **kwargs) -> Expr:
    kwargs.setdefault("evaluate", False)
    return sym_sqrt(expr, **kwargs)


def mul(*exprs: Expr, **kwargs) -> Expr:
    kwargs.setdefault("evaluate", False)
    return Mul(*exprs, **kwargs)
