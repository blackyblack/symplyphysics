from pytest import raises
from sympy import Eq, evaluate
from symplyphysics import Symbol
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import VectorSymbol, VectorNorm as norm
from symplyphysics.core.experimental.solvers import apply, solve_for_vector, vector_equals


def test_apply() -> None:
    x = Symbol("x", real=True)
    y = Symbol("y", real=True)
    eqn = Eq(y, x**2 - x)
    applied_eqn = apply(eqn, lambda side: side + x)
    assert expr_equals(applied_eqn.lhs, x + y)
    assert expr_equals(applied_eqn.rhs, x**2)

    eqn = y - x**2 + x
    applied_eqn = apply(eqn, lambda side: side - x)
    assert expr_equals(applied_eqn.lhs, y - x**2)
    assert expr_equals(applied_eqn.rhs, -1 * x)

    a = VectorSymbol("a")
    b = VectorSymbol("b")
    eqn = Eq(a / norm(a), b / norm(b))
    applied_eqn = apply(eqn, norm)
    assert expr_equals(applied_eqn.lhs, 1)
    assert expr_equals(applied_eqn.rhs, 1)

    eqn = a - b * x
    applied_eqn = apply(eqn, lambda side: side + a)
    assert vector_equals(applied_eqn.lhs, a * 2 - b * x)
    assert vector_equals(applied_eqn.rhs, a)


def test_vector_equals() -> None:
    a = VectorSymbol("a")
    b = VectorSymbol("b")
    c = VectorSymbol("c")

    with evaluate(False):
        lhs = a + b
        rhs = b + a
    assert vector_equals(lhs, rhs)

    with evaluate(False):
        lhs = a + c
        rhs = c + b + a - b
    assert vector_equals(lhs, rhs)

    with evaluate(False):
        lhs = a + b * 3 - a * 2 + a - b - b - b
    assert vector_equals(lhs, 0)


def test_express_atomic() -> None:
    a = VectorSymbol("a")
    b = VectorSymbol("b")
    c = VectorSymbol("c")

    old_expr = a * 2 - b + c
    new_eqn = solve_for_vector(old_expr, b)
    assert vector_equals(new_eqn.lhs, b)
    assert vector_equals(new_eqn.rhs, a * 2 + c)

    new_eqn = solve_for_vector(old_expr, a, reduce_factor=False)
    assert vector_equals(new_eqn.lhs, a * (-2))
    assert vector_equals(new_eqn.rhs, c - b)

    old_eqn = Eq(a * 2, c - b)
    new_eqn = solve_for_vector(old_eqn, b)
    assert vector_equals(new_eqn.lhs, b)
    assert vector_equals(new_eqn.rhs, c - a * 2)

    new_eqn = solve_for_vector(old_eqn, c, reduce_factor=False)
    assert vector_equals(new_eqn.lhs, c)
    assert vector_equals(new_eqn.rhs, a * 2 + b)

    # The expression is not a VectorExpr
    with raises(TypeError):
        solve_for_vector(norm(a), a)

    # The expression does not contain the symbol
    old_expr = a + c
    with raises(ValueError):
        solve_for_vector(old_expr, b)

    # The equation does not contain the symbol
    old_eqn = Eq(a, c)
    with raises(ValueError):
        solve_for_vector(old_eqn, b)
