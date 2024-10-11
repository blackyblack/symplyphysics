from typing import Any
from sympy import Expr, log as sym_log, E, sqrt as sym_sqrt
from sympy.core.parameters import global_parameters

_old_evaluation: bool = True

# Prevents auto processing of expressions. Helps making documentation cleaner.
def disable_sympy_evaluation() -> None:
    _old_evaluation = global_parameters.evaluate
    global_parameters.evaluate = False

def enable_sympy_evaluation() -> None:
    _old_evaluation = global_parameters.evaluate
    global_parameters.evaluate = True

# Restores auto processing of expressions. Allows SymPy work with expressions properly.
def reset_sympy_evaluation() -> None:
    global_parameters.evaluate = _old_evaluation

def log(expr: Expr, base: Expr = E, **kwargs: Any) -> Expr:
    kwargs.setdefault("evaluate", False)
    return sym_log(expr, base, **kwargs)

def sqrt(expr: Expr, **kwargs: Any) -> Expr:
    kwargs.setdefault("evaluate", False)
    return sym_sqrt(expr, **kwargs)
