from sympy import Symbol, Point
from symplyphysics.core.geometry.line import two_point_function
from symplyphysics.core.expr_comparisons import expr_equals

def test_zero_intercept() -> None:
    p1 = Point(0, 0)
    p2 = Point(1, 2)
    x = Symbol("x")
    y = two_point_function(p1, p2, x)
    assert expr_equals(y, 2 * x)


def test_nonzero_intercept() -> None:
    p1 = Point(1, 1)
    p2 = Point(2, 3)
    x = Symbol("x")
    y = two_point_function(p1, p2, x)
    assert expr_equals(y, 2 * x - 1)


def test_parallel_to_y_axis() -> None:
    p1 = Point(1, 1)
    p2 = Point(1, 2)
    x = Symbol("x")
    assert two_point_function(p1, p2, x) is NotImplemented
