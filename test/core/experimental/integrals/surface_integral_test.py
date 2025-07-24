from pytest import raises

from sympy import evaluate, sin, cos, pi, sqrt

from symplyphysics import Symbol, assert_equal, symbols, units
from symplyphysics.core.symbols.symbols import BasicSymbol

from symplyphysics.core.experimental.vectors import VectorDot, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import (AppliedPoint, CARTESIAN,
    CYLINDRICAL, SPHERICAL, CoordinateVector)
from symplyphysics.core.experimental.coordinate_systems.surface import Surface
from symplyphysics.core.experimental.integrals.surface_integral import (SurfaceIntegral,
    INFINITESIMAL_VECTOR_AREA as ds)


def test_non_evaluate() -> None:
    surface_symbol = BasicSymbol("S")

    result = SurfaceIntegral(ds, surface_symbol)
    assert isinstance(result, SurfaceIntegral)
    assert result.args == (ds, surface_symbol)
    assert result.doit() == result

    t1 = Symbol("t1", real=True)
    t2 = Symbol("t2", real=True)

    surface = Surface((t1, t2), AppliedPoint([t1, t2, 0], CARTESIAN))
    assert surface == surface.doit()

    result = SurfaceIntegral(VectorNorm(ds), surface, ((0, 2), (0, 1)), evaluate=False)
    assert isinstance(result, SurfaceIntegral)
    assert result.args == (VectorNorm(ds), surface, ((0, 2), (0, 1)))
    assert result.doit().doit() == 2

    # Undefined bounds

    result = SurfaceIntegral(VectorNorm(ds), surface, evaluate=False)
    assert isinstance(result, SurfaceIntegral)
    assert result.args == (VectorNorm(ds), surface)
    with raises(ValueError):
        result.doit()


def test_cartesian_vector_surface_integral() -> None:
    m = units.meter

    t1 = Symbol("t1", real=True)
    t2 = Symbol("t2", real=True)

    # Surface area of a rectangular surface in the x-y plane

    surface = Surface((t1, t2), AppliedPoint([t1 * m, t2 * m, 0], CARTESIAN))
    result = SurfaceIntegral(VectorNorm(ds), surface, ((0, 2), (-3, 0)))
    assert_equal(result.doit(), 6 * m**2)

    # More complex surface and integrand

    parametrization = AppliedPoint([
        m * sin(t1) * cos(t2),
        m * sin(t1) * sin(t2),
        m * cos(t1),
    ], CARTESIAN)
    surface = Surface((t1, t2), parametrization)

    x, y, z = CARTESIAN.base_scalars
    f = CoordinateVector([x, y, z], CARTESIAN)
    result = SurfaceIntegral(VectorDot(f, ds), surface, ((0, pi / 2), (0, 2 * pi)))
    assert_equal(result.doit(), 2 * pi * m**3)

    result = SurfaceIntegral(VectorDot(f, ds), surface, ((0, pi / 3), (0, 2 * pi)))
    assert_equal(result.doit(), pi * m**3)

    # Undefined bounds

    with raises(ValueError):
        result = SurfaceIntegral(VectorDot(f, ds), surface)


# This is a test for weird `evaluate` behavior. After using once `evaluate(False)`, the
# `vector_diff` function starts returning vectors with new coordinate system instead of
# the original one.
def test_bad_evaluate() -> None:
    with evaluate(False):
        _expression = 2 * symbols.length
    m1 = units.meter

    _, phi1, _ = CYLINDRICAL.base_scalars

    t11 = Symbol("t11", real=True)

    # Surface area of a cylindrical surface, r=1*meter

    surface1 = Surface((phi1, t11), AppliedPoint([m1, phi1, t11 * m1], CYLINDRICAL))
    result1 = SurfaceIntegral(VectorNorm(ds), surface1, ((0, 2 * pi), (0, 1)))
    assert_equal(result1.doit(), 2 * pi * m1**2)


def test_cylindrical_vector_surface_integral() -> None:
    m = units.meter

    rho, phi, z = CYLINDRICAL.base_scalars

    t = Symbol("t", real=True)

    # Surface area of a cylindrical surface, r=1*meter

    surface = Surface((phi, t), AppliedPoint([m, phi, t * m], CYLINDRICAL))
    result = SurfaceIntegral(VectorNorm(ds), surface, ((0, 2 * pi), (0, 1)))
    assert_equal(result.doit(), 2 * pi * m**2)

    # More complex surface and integrand

    parametrization = AppliedPoint([m * sin(t), phi, m * cos(t)], CYLINDRICAL)
    surface = Surface((t, phi), parametrization)

    f = CoordinateVector([rho, 0, z], CYLINDRICAL, surface.parametrization)

    result = SurfaceIntegral(VectorDot(f, ds), surface, ((0, pi / 2), (0, 2 * pi)))
    assert_equal(result.doit(), 2 * pi * m**3)

    result = SurfaceIntegral(VectorDot(f, ds), surface, ((0, pi / 3), (0, 2 * pi)))
    assert_equal(result.doit(), pi * m**3)

    # Undefined bounds

    with raises(ValueError):
        result = SurfaceIntegral(VectorDot(f, ds), surface)


def test_spherical_vector_surface_integral() -> None:
    m = units.meter

    r, theta, phi = SPHERICAL.base_scalars
    t = Symbol("t", real=True)

    # Surface area of a conical surface, theta=pi/4

    surface = Surface((t, phi), AppliedPoint([m * t, pi / 4, phi], SPHERICAL))
    result = SurfaceIntegral(VectorNorm(ds), surface, ((0, 1), (0, pi)))
    assert_equal(result.doit(), sqrt(2) / 4 * pi * m**2)

    # More complex surface and integrand

    surface = Surface((theta, phi), AppliedPoint([m, theta, phi], SPHERICAL))

    f = CoordinateVector([r, 0, 0], SPHERICAL, surface.parametrization)
    result = SurfaceIntegral(VectorDot(f, ds), surface, ((0, pi / 2), (0, 2 * pi)))
    assert_equal(result.doit(), 2 * pi * m**3)

    result = SurfaceIntegral(VectorDot(f, ds), surface, ((0, pi / 3), (0, 2 * pi)))
    assert_equal(result.doit(), pi * m**3)

    # Undefined bounds

    with raises(ValueError):
        result = SurfaceIntegral(VectorDot(f, ds), surface)
