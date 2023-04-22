from sympy import Expr, simplify
from sympy.vector import Vector


def expr_equals(lhs: Expr, rhs: Expr) -> bool:
    val = simplify(lhs - rhs)
    if val == 0:
        return True
    if val == Vector.zero:
        return True
    return False


## SymPy does not allow to compare Abs with non-Abs values so we apply abs() to both sides.
def expr_equals_abs(lhs: Expr, rhs: Expr) -> bool:
    return expr_equals(abs(lhs), abs(rhs))
