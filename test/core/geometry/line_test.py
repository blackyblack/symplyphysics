from sympy import Symbol, Point2D
from pytest import raises
from symplyphysics.core.geometry.line import two_point_function
from symplyphysics.core.expr_comparisons import expr_equals


def test_zero_intercept() -> None:
    p1 = Point2D(0, 0)
    p2 = Point2D(1, 2)
    x = Symbol("x")
    y = two_point_function(p1, p2, x)
    assert expr_equals(y, 2 * x)


def test_nonzero_intercept() -> None:
    p1 = Point2D(1, 1)
    p2 = Point2D(2, 3)
    x = Symbol("x")
    y = two_point_function(p1, p2, x)
    assert expr_equals(y, 2 * x - 1)


def test_parallel_to_y_axis() -> None:
    p1 = Point2D(1, 1)
    p2 = Point2D(1, 2)
    x = Symbol("x")
    with raises(ValueError):
        two_point_function(p1, p2, x)


def test_points_not_unique() -> None:
    pt = Point2D(1, 1)
    x = Symbol("x")
    with raises(ValueError):
        two_point_function(pt, pt, x)
