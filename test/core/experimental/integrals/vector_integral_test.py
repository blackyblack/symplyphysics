from pytest import raises

from sympy import sin, cos, pi

from symplyphysics import Symbol, symbols
from symplyphysics.core.symbols.symbols import BasicSymbol

from symplyphysics.core.experimental.vectors import (clone_as_vector_symbol, vector_diff,
    clone_as_vector_function)
from symplyphysics.core.experimental.coordinate_systems import (CARTESIAN, CoordinateVector,
    QuantityCoordinateVector)
from symplyphysics.core.experimental.integrals.vector_integral import integrate_coordinate_vectors
from symplyphysics.core.experimental.solvers import vector_equals


def test_integrate_single_vector() -> None:
    t = Symbol("t", real=True)
    s = Symbol("s", real=True)

    expr = CoordinateVector([sin(t), 0, cos(t)], CARTESIAN)
    assert isinstance(expr, CoordinateVector)

    result = integrate_coordinate_vectors(expr, t)
    assert isinstance(result, CoordinateVector)
    assert vector_equals(result, CoordinateVector([-cos(t), 0, sin(t)], CARTESIAN))

    result = integrate_coordinate_vectors(expr, t, t)
    assert isinstance(result, CoordinateVector)
    assert vector_equals(result, CoordinateVector([-sin(t), 0, -cos(t)], CARTESIAN))

    result = integrate_coordinate_vectors(expr, s)
    assert vector_equals(result, CoordinateVector.from_expr(s * expr))

    result = integrate_coordinate_vectors(expr, (t, 0, pi))
    assert isinstance(result, CoordinateVector)
    assert vector_equals(result, CoordinateVector([2, 0, 0], CARTESIAN))

    result = integrate_coordinate_vectors(expr, (t, 0, 2 * pi))
    assert result == 0


def test_integrate_quantity_coordinate_vector() -> None:
    t = Symbol("t", real=True)
    expr = QuantityCoordinateVector([4, 0, 0], CARTESIAN)
    result = integrate_coordinate_vectors(expr, t)
    assert vector_equals(result, CoordinateVector([4 * t, 0, 0], CARTESIAN))

    result = integrate_coordinate_vectors(expr, (t, -1, 1))
    assert isinstance(result, QuantityCoordinateVector)
    assert vector_equals(result, QuantityCoordinateVector([8, 0, 0], CARTESIAN))


def test_integrate_vector_sum() -> None:
    t = Symbol("t", real=True)

    p = BasicSymbol("P")
    q = BasicSymbol("Q")

    u = CoordinateVector([t, 0, 0], CARTESIAN, p)
    v = CoordinateVector([0, t, 0], CARTESIAN, q)
    expr = u + v

    result = integrate_coordinate_vectors(expr, t)
    expected = (CoordinateVector([t**2 / 2, 0, 0], CARTESIAN, p) +
        CoordinateVector([0, t**2 / 2, 0], CARTESIAN, q))
    assert vector_equals(result, expected)


def test_bad_input() -> None:
    t = symbols.time
    r = clone_as_vector_symbol(symbols.distance_to_origin)
    v = clone_as_vector_function(symbols.speed, (t,))

    with raises(TypeError):
        integrate_coordinate_vectors(r, t)
    with raises(TypeError):
        integrate_coordinate_vectors(v(t), t)
    with raises(TypeError):
        integrate_coordinate_vectors(vector_diff(v(t), t), t)

    expr = CoordinateVector([1, 0, 0], CARTESIAN)

    with raises(ValueError):
        integrate_coordinate_vectors(expr, r)

    with raises(ValueError):
        integrate_coordinate_vectors(expr, t, r)
