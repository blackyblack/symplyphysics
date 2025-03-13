from sympy import solve as sym_solve, Eq, Expr

from ..vectors import VectorExpr, VectorNorm, VectorAdd, norm

accepted_types = (VectorNorm,)


def solve_for_expr(f: Eq | Expr | VectorExpr, x: Expr) -> list[Expr]:
    f = f.doit()

    if isinstance(f, Eq):
        f = f.lhs - f.rhs

    if isinstance(f, Expr):
        return sym_solve(f, x)

    return _solve_vector_for_expr(f, x)


def _solve_vector_for_expr(f: VectorExpr, x: Expr) -> list[Expr]:
    if isinstance(x, VectorNorm):
        return _solve_vector_for_norm(f, x)

    # Assume that `x` is not a scalar-valued vector function (such as `norm`).

    combinations = f.as_symbol_combination()

    # Count how many factors contain `x`
    indeces: list[int] = []
    for i, (_, s) in enumerate(combinations):
        if s.has(x):
            indeces.append(i)

    n = len(indeces)

    if n == 0:
        return []

    # (2) `v * (-k(x)) = w` where `w` does not depend on `x`
    if n == 1:
        v, k = combinations[i][0], combinations[i][1]
        ws = combinations[:i] + combinations[i + 1:]
        n_ws = len(ws)

        # (2.1) `v * k(x) = 0` => solve `k(x) = 0` for `x` (another option is `v = 0`)
        if n_ws == 0:
            return sym_solve(k, x)

        w = VectorAdd(*((u * s) for u, s in ws))

        # (2.2) `v * (-k(x)) = w` => solve `norm(v) * abs(k(x)) = norm(w)` for `x`
        equation = Eq(norm(v) * abs(k), norm(w))

        solutions = sym_solve(equation, x)
        return solutions

    # (3) `n >= 2`
    raise NotImplementedError("Dot product has not been implemented yet for this.")


def _solve_vector_for_norm(f: VectorExpr, x: VectorNorm) -> list[Expr]:
    raise NotImplementedError
