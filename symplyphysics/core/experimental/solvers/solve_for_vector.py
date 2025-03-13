from typing import Optional
from sympy import solve as sym_solve, Eq, ask, Q

from ..vectors import VectorExpr, VectorSymbol, ZERO, VectorAdd, VectorNorm, norm


def solve_for_vector_symbol(f: Eq | VectorExpr, x: VectorSymbol) -> list[VectorExpr]:
    if isinstance(f, Eq):
        f = f.lhs - f.rhs

    if not isinstance(f, VectorExpr):
        raise TypeError(f"Solving for a vector requires a vector expression, got {f}.")

    # (1) `combination = x * k + sum(v_i * k_i, (i, 1, n))` where `k` might be `0`.
    combination = f.doit().as_symbol_combination()

    # By the construction of `.as_symbol_combination`, `x` can appear in at most one LHS within
    # `combination`.
    index: Optional[int] = None

    for i, (v, _) in enumerate(combination):
        if v == x:
            index = i
            break

    # `k = 0`, see (1).
    if index is None:
        return []

    # (2) `x * factor = sum(v_i * k_i, (i, 1, n))` where `factor = -k` and `rhs[i] = v_i * k_i`.
    factor = -combination[i][1]
    rhs = combination[:i] + combination[i + 1:]
    n_rhs = len(rhs)

    # (3.1) `rhs` is empty, therefore `Eq(x, 0) | Eq(k, 0)`.
    if n_rhs == 0:
        return [ZERO]

    # (3.2) `n_rhs >= 2`, check if the vector symbol appears in any of the scalar factors.
    if n_rhs >= 2 and any(s.has(x) for _, s in rhs):
        raise NotImplementedError("Vector is contained in more than two terms.")  # TODO

    # (3.3) `x * factor = y` where `y` may have several terms.
    y = VectorAdd(*(v * s for v, s in rhs))

    # (3.3.1) `factor` does not depend on `x`, therefore `x = y / factor`.
    if not factor.has(x):
        return [y / factor]

    # (3.3.2) `factor` depends on `norm(x)`, i.e. `factor(norm(x)) * x = y`. In this case, we can
    # solve the equation `abs(factor(norm(x))) * norm(x) * norm(y)` for `norm(x)`. As a result, we
    # can write `x = y / factor(norm(x))` where `norm(x)` has been substituted for the
    # aforementioned solution(s).
    if factor.has(VectorNorm(x)):
        equation = Eq(abs(factor) * VectorNorm(x), norm(y))

        solved_norms = []
        for e in sym_solve(equation, VectorNorm(x)):
            if ask(Q.nonnegative(e)):
                solved_norms.append(e)

        solved_factors = [factor.subs(VectorNorm(x), e) for e in solved_norms]
        return [y / e for e in solved_factors]

    raise NotImplementedError("Unknown vector classes in input.")
