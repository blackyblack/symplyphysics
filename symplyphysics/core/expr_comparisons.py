from sympy import Expr, simplify

def expr_equals(lhs: Expr, rhs: Expr) -> bool:
    return simplify(lhs - rhs) == 0

## SymPy does not allow to compare Abs with non-Abs values so we apply abs() to both sides.
def expr_equals_abs(lhs: Expr, rhs: Expr) -> bool:
    return expr_equals(abs(lhs), abs(rhs))
