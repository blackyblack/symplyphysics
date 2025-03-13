from pytest import raises
from sympy import Eq, Symbol as SymSymbol
from symplyphysics import symbols, units, clone_as_symbol
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import VectorSymbol, norm
from symplyphysics.core.experimental.solvers.solve_for_expr import solve_for_expr


def test_solve_for_expr() -> None:
    # One-term expression

    v = VectorSymbol("v")
    k = SymSymbol("k")
    eqn = v * k * (k - 1)

    solved_ks = sorted(solve_for_expr(eqn, k))
    assert expr_equals(solved_ks[0], 0)
    assert expr_equals(solved_ks[1], 1)

    # Two-term expression.

    f = VectorSymbol("F", units.force)
    m = clone_as_symbol(symbols.mass, positive=True)
    a = VectorSymbol("a", units.acceleration)
    eqn = Eq(f, a * m)

    solved_m, = solve_for_expr(eqn, m)
    assert expr_equals(solved_m, norm(f) / norm(a))

    with raises(NotImplementedError):
        solve_for_expr(eqn, norm(f))

    with raises(NotImplementedError):
        solve_for_expr(eqn, norm(a))

    # Expression of interest is in several terms.

    r = VectorSymbol("r", units.length)
    r0 = VectorSymbol("r_0", units.length)
    v0 = VectorSymbol("v_0", units.velocity)
    t = clone_as_symbol(symbols.time, positive=True)
    eqn = Eq(r, r0 + v0 * t + a * t**2 / 2)

    with raises(NotImplementedError):
        solve_for_expr(eqn, t)
