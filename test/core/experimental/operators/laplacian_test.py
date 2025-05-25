from pytest import raises

from sympy import Symbol as SymSymbol, sin, cos, tan

from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import VectorSymbol
from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    CylindricalCoordinateSystem, SphericalCoordinateSystem)
from symplyphysics.core.experimental.coordinate_systems.vector import CoordinateVector
from symplyphysics.core.experimental.coordinate_systems.scalar import CoordinateScalar
from symplyphysics.core.experimental.operators import VectorLaplacian


def test_non_coordinate_scalar() -> None:
    assert expr_equals(VectorLaplacian(0), 0)
    assert expr_equals(VectorLaplacian(1), 0)

    x = SymSymbol("x", real=True)
    assert expr_equals(VectorLaplacian(x), 0)

    v = VectorSymbol("v")
    with raises(ValueError):
        _ = VectorLaplacian(v)

    cartesian = CartesianCoordinateSystem()
    v = CoordinateVector([1, 2, 3], cartesian)
    with raises(ValueError):
        _ = VectorLaplacian(v)


def test_non_evaluate() -> None:
    assert expr_equals(VectorLaplacian(0, evaluate=False), 0)

    result = VectorLaplacian(1, evaluate=False)
    assert isinstance(result, VectorLaplacian)
    assert result.args == (1,)
    assert expr_equals(result, 0)  # evaluate

    x = SymSymbol("x", real=True)
    result = VectorLaplacian(2 * x, evaluate=False)
    assert isinstance(result, VectorLaplacian)
    assert result.args == (2 * x,)
    assert expr_equals(result, 0)  # evaluate

    v = VectorSymbol("v")
    with raises(ValueError):
        _ = VectorLaplacian(v, evaluate=False)

    cartesian = CartesianCoordinateSystem()
    v = CoordinateVector([1, 2, 3], cartesian)
    with raises(ValueError):
        _ = VectorLaplacian(v, evaluate=False)


def test_cartesian_system() -> None:
    cartesian = CartesianCoordinateSystem()

    x, y, z = cartesian.base_scalars

    scalar = CoordinateScalar(sin(x + y + z), cartesian)

    result = VectorLaplacian(scalar)
    assert isinstance(result, CoordinateScalar)
    assert expr_equals(result.scalar, -3 * sin(x + y + z))
    assert result.system == cartesian
    assert result.point == scalar.point

    scalar = CoordinateScalar(x * y * z, cartesian)
    assert expr_equals(VectorLaplacian(scalar), 0)


def test_cylindrical_system() -> None:
    cylindrical = CylindricalCoordinateSystem()

    rho, phi, z = cylindrical.base_scalars
    p = SymSymbol("P")

    scalar: CoordinateScalar = CoordinateScalar(rho**2 + z * sin(phi), cylindrical, p)
    result = VectorLaplacian(scalar)
    assert isinstance(result, CoordinateScalar)
    assert expr_equals(result.scalar, 4 - z * sin(phi) / rho**2)
    assert result.system == scalar.system
    assert result.point == scalar.point

    scalar = CoordinateScalar((rho + z)**2 * cos(phi) / 3, cylindrical, p)
    result = VectorLaplacian(scalar)
    assert isinstance(result, CoordinateScalar)
    assert expr_equals(result.scalar.simplify(), (5 - (z / rho)**2) * (cos(phi) / 3))
    assert result.system == scalar.system
    assert result.point == scalar.point


def test_spherical_system() -> None:
    spherical = SphericalCoordinateSystem()

    r, theta, phi = spherical.base_scalars
    p = SymSymbol("P")

    scalar: CoordinateScalar = CoordinateScalar(r * sin(theta) * cos(phi), spherical, p)
    result = VectorLaplacian(scalar)
    assert isinstance(result, CoordinateScalar)
    assert expr_equals(
        result.scalar,
        0,
    )

    scalar = CoordinateScalar(r**2 * cos(theta + phi), spherical, p)
    result = VectorLaplacian(scalar)
    assert isinstance(result.simplify(), CoordinateScalar)
    assert expr_equals(
        result.scalar,
        5 * cos(theta + phi) - sin(theta + phi) / tan(theta) - cos(theta + phi) / sin(theta)**2,
    )
