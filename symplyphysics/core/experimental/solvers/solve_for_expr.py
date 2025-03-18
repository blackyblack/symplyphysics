from typing import Any
from sympy import solve as sym_solve, Eq, Expr, Symbol

from ..vectors import VectorExpr, VectorNorm, VectorAdd, norm, VectorSymbol


class FakeVectorSymbol(Symbol):
    pass


def solve_for_expr(f: Eq | Expr | VectorExpr, x: Expr, **kwargs: Any) -> list[Expr]:
    if isinstance(f, Eq):
        f = f.lhs - f.rhs

    if isinstance(f, Expr):
        return sym_solve(f, x, **kwargs)

    if not isinstance(f, VectorExpr):
        raise TypeError(f"Unknown type {type(f).__name__} for an equation.")

    if isinstance(x, VectorNorm):
        return solve_for_norm(f, x)

    raise NotImplementedError


def solve_for_norm(f: VectorExpr, x: VectorNorm) -> list[Expr]:
    w = x.argument
    if not isinstance(w, VectorSymbol):
        raise ValueError("...")  # TODO

    combination = f.as_symbol_combination()

    i_w = None
    i_ss = []

    for j, (v, s) in enumerate(combination):
        if v == w:
            i_w = j

        if s.has(w):
            i_ss.append(j)

    def get_rest(*excluded: int) -> VectorExpr:
        excluded = set(excluded)
        pairs = (p for (i, p) in enumerate(combination) if i not in excluded)
        return VectorAdd(*(v * s for v, s in pairs))

    if i_w is None:
        # `0 = sum(v_i * s_i(norm(w)) | 1 <= i <= n)` where for all `v_i`, `v_i != w`

        if not i_ss:
            # none of `s_i` depend on `norm(w)`
            return []

        if len(i_ss) > 2:
            # Attempt to convert to the form `a = b * f(norm(w))` where `a, b` are vectors independent of `norm(w)`
            raise NotImplementedError("...")  # TODO

        # `0 = v * s(norm(w)) + rest`
        # `=> norm(v)^2 * s(norm(w))^2 = norm(rest)^2`
        i_s = i_ss[0]
        rest = get_rest(i_s)
        v, s = combination[i_s]
        eqn = Eq(norm(v)**2 * abs(s)**2, norm(rest))
        return sym_solve(eqn, x)

    if not i_ss:
        # `0 = w * s + rest` where neither `s` nor `rest` depend on `norm(w)`
        # `=> w * s = -rest`
        # `=> norm(w) * abs(s) = norm(rest)`
        # `=> norm(w) = norm(rest) / abs(s)`

        rest = get_rest(i_w)
        return norm(rest) / abs(combination[i_w][1])

    if i_w in i_ss:
        pass

    # Attempt to convert to the from `w * s = v * f(norm(w))` where scalar `s` and vector `v` don't depend on `norm(w)`
    raise NotImplementedError("...")  # TODO
