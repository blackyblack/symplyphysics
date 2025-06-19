from pytest import raises

from sympy import S, log, sin, cos, pi, sqrt, Symbol as SymSymbol

from symplyphysics import Symbol, assert_equal, units
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.vectors import VectorDot
from symplyphysics.core.experimental.coordinate_systems import (AppliedPoint, CARTESIAN,
    CoordinateVector)
from symplyphysics.core.experimental.coordinate_systems.curve import Curve
from symplyphysics.core.experimental.integrals.line_integral import (LineIntegral,
    INFINITESIMAL_ARC_LENGTH as ds, INFINITESIMAL_DISPLACEMENT as dr)

x, y, z = CARTESIAN.base_scalars


def test_non_evaluate() -> None:
    curve_symbol = SymSymbol("C")

    result = LineIntegral(ds, curve_symbol)
    assert isinstance(result, LineIntegral)
    assert result.args == (ds, curve_symbol)
    assert result.doit() == result

    t = Symbol("t", positive=True)
    curve = Curve(t, AppliedPoint([t, 0, 0], CARTESIAN))
    assert curve == curve.doit()

    result = LineIntegral(ds, curve, (0, 1), evaluate=False)
    assert isinstance(result, LineIntegral)
    assert result.args == (ds, curve, (0, 1))
    assert result.doit().doit() == 1

    # Undefined bounds
    result = LineIntegral(ds, curve, evaluate=False)
    assert isinstance(result, LineIntegral)
    assert result.args == (ds, curve)
    with raises(ValueError):
        result.doit()


def test_cartesian_scalar_line_integral() -> None:
    t = Symbol("t", positive=True)

    # The curve is just the x-axis

    curve = Curve(t, AppliedPoint([t, 0, 0], CARTESIAN))

    result = LineIntegral(ds, curve, (0, 1))
    assert expr_equals(result, 1)

    result = LineIntegral(x * ds, curve, (0, 1))
    assert expr_equals(result, S.Half)

    result = LineIntegral((1 + x**2) * ds, curve, (0, 1))
    assert expr_equals(result, 1 + S(1) / 3)

    # A more intricate curve

    curve = Curve(t, AppliedPoint([t**2, 1 - t, log(t)], CARTESIAN))
    result = LineIntegral(ds / sqrt(x), curve, (1, 2))
    assert_equal(result.doit(), 2.19, relative_tolerance=2e-3)

    # Undefined bounds
    with raises(ValueError):
        LineIntegral(ds / sqrt(x), curve)


def test_cartesian_vector_line_integral() -> None:
    t = Symbol("t", real=True)

    curve = Curve(t, AppliedPoint([t, t, t], CARTESIAN))

    field = CoordinateVector([
        sin(pi * y / 2),
        cos(pi * x / 2),
        0,
    ], CARTESIAN)

    result = LineIntegral(VectorDot(field, dr), curve, (-1, 2))
    assert expr_equals(result, 4 / pi)

    curve = Curve(t, AppliedPoint([0, sin(t), cos(t)], CARTESIAN))

    field = CoordinateVector([t, -t**2, -t], CARTESIAN)
    result = LineIntegral(VectorDot(field, dr), curve, (0, 2 * pi))
    assert expr_equals(result, -6 * pi)

    # Undefined bounds
    with raises(ValueError):
        LineIntegral(VectorDot(field, dr), curve)


def test_cartesian_vector_line_integral_with_quantities() -> None:
    t = Symbol("t", real=True)

    field = CoordinateVector([y, 0, x + z], CARTESIAN)
    curve = Curve(t, AppliedPoint([units.meter * cos(t), units.meter * sin(t), 0], CARTESIAN))
    result = LineIntegral(VectorDot(field, dr), curve, (0, pi / 2))
    assert_equal(result.doit(), -pi / 4 * units.meter**2)

    field = CoordinateVector([0, -1 * units.newton * (units.meter / y)**2, 0], CARTESIAN)

    curve = Curve(t, AppliedPoint([units.meter * t, units.meter * t, 0], CARTESIAN))
    result = LineIntegral(VectorDot(field, dr), curve, (1, 2))
    assert_equal(result.doit(), -0.5 * units.joule)

    # y = 5
    curve = Curve(t, AppliedPoint([units.meter * t, units.meter * 5, 0], CARTESIAN))
    result = LineIntegral(VectorDot(field, dr), curve, (1, 2))
    assert_equal(result.doit(), 0)

    # x = 5
    curve = Curve(t, AppliedPoint([units.meter * 5, units.meter * t, 0], CARTESIAN))
    result = LineIntegral(VectorDot(field, dr), curve, (1, 2))
    assert_equal(result.doit(), -0.5 * units.joule)

    # Change direction to top-down
    result = LineIntegral(VectorDot(field, dr), curve, (2, 1))
    assert_equal(result.doit(), 0.5 * units.joule)


# TODO: add cylindrical and spherical tests
