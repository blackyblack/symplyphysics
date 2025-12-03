from typing import Callable
from sympy import Basic, S, Eq, sympify, Add, Expr

from ..vectors import (
    VectorSymbol,
    VectorNorm as norm,
    is_vector_expr,
    into_terms,
    split_factor,
)


def apply(eqn: Basic, f: Callable[[Basic], Basic]) -> Eq:
    """
    Applies `f` to both sides of the equation `eqn`, if it is an equality, else treats `eqn` as the
    left-hand side and `zero` as the right-hand side and applies `f` to them. Returns an equation
    object.
    """

    eqn = sympify(eqn, strict=True)

    if isinstance(eqn, Eq):
        lhs = eqn.lhs
        rhs = eqn.rhs
    else:
        lhs = eqn
        rhs = S.Zero

    return Eq(f(lhs), f(rhs), evaluate=False)


def vector_equals(lhs: Expr, rhs: Expr) -> bool:
    """Checks the equality of two vector expressions."""

    diff = (lhs - rhs).simplify()
    return bool(norm(diff) == 0)


def solve_for_vector(
    expr: Eq | Expr,
    atomic: VectorSymbol,
) -> Expr:
    """
    Viewing `expr` as a linear combination of vectors, expresses the `atomic` vector using the rest
    of the vectors comprising that linear combination.

    To elaborate, suppose `expr` is written in the form `sum(k_i * v_i, i)` where `{k_i}` are
    scalars and `{v_i}` are (atomic) vectors (which can always be done for vector expressions).
    Then, if `atomic` is found under index `j` within `{v_i}`, rewrites `expr` in the form `atomic
    = sum(-k_i * v_i, i ≠ j)`.

    Note that whether any of `k_i` depend on `atomic` or not is not taken into account, so it does
    not truly "solve" w.r.t. `atomic` unless none of `k_i` depend on `atomic`, in which case the
    result of `solve_for_vector` is indeed the solution of `expr` w.r.t. `atomic`. If some of
    `{k_i}` do depend on `atomic`, though, a scalar-reducing operation, such as taking the dot
    product of both sides of the equation with a known vector, might help find a solution.

    Raises:
        TypeError: If `expr` is not a vector expression or vector equation.
        ValueError: If `atomic` does not appear in the linear combination in `expr`.
    """

    if isinstance(expr, Eq):
        expr = expr.lhs - expr.rhs

    if not is_vector_expr(expr):
        raise ValueError(f"Expected '{expr}' to be a vector.")

    combination = tuple(split_factor(term) for term in into_terms(expr))
    i = None

    for j, (v, _) in enumerate(combination):
        if vector_equals(v, atomic):
            i = j
            break

    if i is None:
        raise ValueError(f"The expression {expr} does not contain the symbol {atomic}.")

    combination_rhs = combination[:i] + combination[i + 1:]
    scale = combination[i][1]

    rhs = Add(*(v * (-1 * s / scale) for v, s in combination_rhs))
    return rhs
