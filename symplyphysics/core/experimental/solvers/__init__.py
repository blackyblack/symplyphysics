from typing import Callable, Any
from sympy import Basic, S, Eq, sympify, solve as sym_solve
from ..vectors import VectorExpr, VectorSymbol, norm, VectorAdd


def apply(eqn: Any, f: Callable[[Basic], Basic], *, zero: Basic = S.Zero) -> Eq:
    """
    Applies `f` to both sides of the equation `eqn`, if it is an equality, else treats `f` as the
    left-hand side and `zero` as the right-hand side and applies `f` to them. Returns an equation
    object.

    Note that if `eqn` is neither `Eq` nor `Expr`, the `zero` should be specified to avoid errors.
    """

    eqn = sympify(eqn, strict=True)

    if isinstance(eqn, Eq):
        lhs = eqn.lhs
        rhs = eqn.rhs
    else:
        lhs = eqn
        rhs = zero

    return Eq(f(lhs), f(rhs), evaluate=False)


def solve_into_eq(f: Any, symbol: Any, **flags: Any) -> list[Eq]:
    flags["dict"] = True
    solution = sym_solve(f, symbol, **flags)[0]
    return [Eq(lhs, rhs) for lhs, rhs in solution.items()]


def vector_equals(lhs: VectorExpr, rhs: VectorExpr) -> bool:
    """Checks the equality of two vector expressions."""

    diff = lhs - rhs
    return bool(norm(diff) == 0)


def express_atomic(
    expr: Eq | VectorExpr,
    atomic: VectorSymbol,
    reduce_factor: bool = True,
) -> Eq:
    """
    Rewrites `expr` so that `atomic` appears in the LHS. Note that only the symbols that make up
    the linear combination in `expr` are taken into account when searching for `atomic`.

    Args:
        expr: A vector expression or vector equation.
        atomic: The vector symbol to be expressed.
        reduce_factor:
            When set to `True`, reduces both sides of the equation by the factor to which `atomic`
            is multiplied.

    Returns:
        An equation with `atomic` in the LHS, or `None` if `atomic` is not present in `expr`.

    Raises:
        TypeError: If `expr` is not a vector expression or vector equation.
        ValueError: If `atomic` does not appear in the linear combination in `expr`.
    """

    if isinstance(expr, Eq):
        expr = expr.lhs - expr.rhs

    if not isinstance(expr, VectorExpr):
        raise TypeError("...")

    combination = expr.as_symbol_combination()
    i = None

    for j, (v, _) in enumerate(combination):
        if vector_equals(v, atomic):
            i = j
            break

    if i is None:
        raise ValueError("...")

    combination_rhs = combination[:i] + combination[i + 1:]
    scale = combination[i][1]

    if reduce_factor:
        rhs = VectorAdd(*(v * (-1 * s / scale) for v, s in combination_rhs))
        return Eq(atomic, rhs)

    rhs = VectorAdd(*(v * s for v, s in combination_rhs))
    return Eq(atomic * (-1 * scale), rhs)
