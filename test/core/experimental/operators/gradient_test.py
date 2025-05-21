from typing import Sequence, Any
from pytest import raises

from sympy import sin, cos, Symbol as SymSymbol, Expr, ImmutableMatrix

from symplyphysics.core.experimental.vectors import VectorSymbol
from symplyphysics.core.experimental.coordinate_systems import (BaseCoordinateSystem,
    CartesianCoordinateSystem, CylindricalCoordinateSystem, SphericalCoordinateSystem)
from symplyphysics.core.experimental.coordinate_systems.vector import CoordinateVector
from symplyphysics.core.experimental.coordinate_systems.scalar import CoordinateScalar
from symplyphysics.core.experimental.operators import VectorGradient
from symplyphysics.core.experimental.solvers import vector_equals


def test_non_coordinate_scalar() -> None:
    assert vector_equals(VectorGradient(0), 0)
    assert vector_equals(VectorGradient(1), 0)

    x = SymSymbol("x", real=True)
    assert vector_equals(VectorGradient(x), 0)  # not a CoordinateScalar

    v = VectorSymbol("v")
    with raises(ValueError):
        _ = VectorGradient(v)

    cartesian = CartesianCoordinateSystem()
    v = CoordinateVector([1, 2, 3], cartesian)
    with raises(ValueError):
        _ = VectorGradient(v)


def test_non_evaluate() -> None:
    assert vector_equals(VectorGradient(0, evaluate=False), 0)

    grad = VectorGradient(1, evaluate=False)
    assert isinstance(grad, VectorGradient)
    assert grad.args == (1,)
    assert vector_equals(grad, 0)  # evaluate

    x = SymSymbol("x", real=True)
    grad = VectorGradient(x, evaluate=False)
    assert isinstance(grad, VectorGradient)
    assert grad.args == (x,)
    assert vector_equals(grad, 0)  # evaluate

    v = VectorSymbol("v")
    with raises(ValueError):
        _ = VectorGradient(v, evaluate=False)

    cartesian = CartesianCoordinateSystem()
    v = CoordinateVector([1, 2, 3], cartesian)
    with raises(ValueError):
        _ = VectorGradient(v, evaluate=False)

    x = cartesian.base_scalars[0]
    s = CoordinateScalar(x, cartesian)
    grad = VectorGradient(s, evaluate=False)
    assert grad.args == (s,)
    assert vector_equals(grad.doit(), CoordinateVector([1, 0, 0], cartesian))


def check(system: BaseCoordinateSystem, value: Expr, components: Sequence[Any]) -> None:
    p = SymSymbol("P")

    scalar = CoordinateScalar(value, system, p)

    gradient = VectorGradient(scalar)

    expected = CoordinateVector(components, system, p)

    assert vector_equals(gradient, expected)


def test_cartesian_system() -> None:
    system = CartesianCoordinateSystem()

    x, y, z = system.base_scalars

    value = sin(x + y + z)
    check(system, value, ImmutableMatrix.ones(3, 1) * cos(x + y + z))

    value = (x**2 + y**2 + z**2) / 2
    check(system, value, [x, y, z])

    value = x * y + y * z + z * x
    check(system, value, [y + z, x + z, x + y])


def test_cylindrical_system() -> None:
    system = CylindricalCoordinateSystem()

    rho, phi, z = system.base_scalars

    value = rho * sin(phi) - z
    check(system, value, [sin(phi), cos(phi), -1])

    value = (rho**2 + z**2) / 2
    check(system, value, [rho, 0, z])


def test_spherical_system() -> None:
    system = SphericalCoordinateSystem()

    r, theta, phi = system.base_scalars

    value = r * sin(theta - phi)
    check(system, value, [
        sin(theta - phi),
        cos(theta - phi),
        -cos(theta - phi) / sin(theta),
    ])

    value = r**2 / 2
    check(system, value, [r, 0, 0])
