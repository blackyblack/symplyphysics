from sympy import symbols, sin, cos, pi
from symplyphysics.core.expr_comparisons import expr_equals, expr_equals_abs


def test_basic_comparison():
    x1, x2 = symbols("x1 x2")
    assert expr_equals(x1, x1)
    assert not expr_equals(x1, x2)


def test_trig_comparison():
    x1 = symbols("x1")
    assert expr_equals(sin(pi / 2 - x1), cos(x1))
    assert not expr_equals(sin(x1), cos(x1))


def test_numeric_comparison():
    x1 = 1
    x2 = 2
    assert expr_equals(x1, 1)
    assert expr_equals(3, x1 + x2)


def test_addition_comparison():
    x1, x2 = symbols("x1 x2")
    assert expr_equals(x1 + x2, x2 + x1)


def test_basic_abs_comparison():
    x1, x2 = symbols("x1 x2")
    assert expr_equals_abs(x1, x1)
    assert expr_equals_abs(-x1, x1)
    assert expr_equals_abs(x1, -x1)
    assert expr_equals_abs(-x1, -x1)
    assert not expr_equals_abs(x1, x2)
