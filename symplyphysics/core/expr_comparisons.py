from typing import SupportsAbs, Any
from sympy import simplify
from sympy.vector import Vector


## Do not try to limit type of the input parameters. Allow any object to
## be compared, if it can.
def expr_equals(lhs: Any, rhs: Any) -> bool:
    val = simplify(lhs - rhs)
    if val == 0:
        return True
    if val == Vector.zero:
        return True
    return False


## SymPy does not allow to compare Abs with non-Abs values so we apply abs() to both sides.
def expr_equals_abs(lhs: SupportsAbs, rhs: SupportsAbs) -> bool:  # type: ignore[type-arg]
    return expr_equals(abs(lhs), abs(rhs))
