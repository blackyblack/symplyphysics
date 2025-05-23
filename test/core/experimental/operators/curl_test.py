from pytest import raises

from sympy import ImmutableMatrix, Symbol as SymSymbol, sin, cos, tan

from symplyphysics.core.experimental.vectors import VectorSymbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    CylindricalCoordinateSystem, SphericalCoordinateSystem)
from symplyphysics.core.experimental.coordinate_systems.vector import CoordinateVector
from symplyphysics.core.experimental.coordinate_systems.scalar import CoordinateScalar
from symplyphysics.core.experimental.operators import VectorCurl
from symplyphysics.core.experimental.solvers import vector_equals


def test_non_coordinate_scalar() -> None:
    assert vector_equals(VectorCurl(0), 0)

    with raises(ValueError):
        assert vector_equals(VectorCurl(1), 0)

    x = SymSymbol("x", real=True)
    with raises(ValueError):
        assert vector_equals(VectorCurl(x), 0)

    v = VectorSymbol("v")
    assert VectorCurl(v) == 0

    w = VectorSymbol("w")
    assert VectorCurl(v + w) == 0
    assert VectorCurl(VectorCross(v, w)) == 0

    cartesian = CartesianCoordinateSystem()
    v = CoordinateScalar(1, cartesian)
    with raises(ValueError):
        _ = VectorCurl(v)


def test_non_evaluate() -> None:
    assert vector_equals(VectorCurl(0, evaluate=False), 0)

    with raises(ValueError):
        assert vector_equals(VectorCurl(1, evaluate=False), 0)

    x = SymSymbol("x", real=True)
    with raises(ValueError):
        assert vector_equals(VectorCurl(x), 0)

    v = VectorSymbol("v")
    curl = VectorCurl(v, evaluate=False)
    assert isinstance(curl, VectorCurl)
    assert curl.args == (v,)
    assert vector_equals(curl, 0)  # evaluate

    w = VectorSymbol("w")
    curl = VectorCurl(v + w, evaluate=False)
    assert isinstance(curl, VectorCurl)
    assert curl.args == (v + w,)
    assert vector_equals(curl, 0)  # evaluate

    curl = VectorCurl(VectorCross(v, w), evaluate=False)
    assert isinstance(curl, VectorCurl)
    assert curl.args == (VectorCross(v, w),)
    assert vector_equals(curl, 0)  # evaluate

    cartesian = CartesianCoordinateSystem()
    vector = CoordinateVector([1, 1, 1], cartesian)
    curl = VectorCurl(vector, evaluate=False)
    assert isinstance(curl, VectorCurl)
    assert curl.args == (vector,)
    assert vector_equals(curl, 0)  # evaluate


def test_cartesian_system() -> None:
    cartesian = CartesianCoordinateSystem()

    x, y, z = cartesian.base_scalars

    vector = CoordinateVector([z, x, y], cartesian)
    assert vector_equals(
        VectorCurl(vector),
        CoordinateVector(ImmutableMatrix.ones(3, 1), cartesian),
    )

    vector = CoordinateVector([x, y, z], cartesian)
    assert vector_equals(VectorCurl(vector), 0)


def test_cylindrical_system() -> None:
    cylindrical = CylindricalCoordinateSystem()

    rho, phi, z = cylindrical.base_scalars
    p = SymSymbol("P")

    vector = CoordinateVector([rho, z * sin(phi), z], cylindrical, p)
    assert vector_equals(
        VectorCurl(vector),
        CoordinateVector([-sin(phi), 0, z * sin(phi) / rho], cylindrical, p),
    )

    vector = CoordinateVector([z * sin(phi), rho, rho * cos(phi)], cylindrical, p)
    expected = CoordinateVector([
        -sin(phi),
        sin(phi) - cos(phi),
        2 - z * cos(phi) / rho,
    ], cylindrical, p)
    assert vector_equals(VectorCurl(vector), expected)


def test_spherical_system() -> None:
    spherical = SphericalCoordinateSystem()

    r, theta, phi = spherical.base_scalars
    p = SymSymbol("P")

    vector = CoordinateVector([r * sin(theta), r * sin(phi), r], spherical, p)
    expected = CoordinateVector([
        (cos(theta) - cos(phi)) / sin(theta),
        -2,
        2 * sin(phi) - cos(theta),
    ], spherical, p)
    assert vector_equals(VectorCurl(vector), expected)

    vector = CoordinateVector([r, 0, 0], spherical, p)
    assert vector_equals(VectorCurl(vector), 0)

    vector = CoordinateVector([0, r, 0], spherical, p)
    assert vector_equals(
        VectorCurl(vector),
        CoordinateVector([0, 0, 2], spherical, p),
    )

    vector = CoordinateVector([0, 0, r], spherical, p)
    assert vector_equals(
        VectorCurl(vector),
        CoordinateVector([1 / tan(theta), -2, 0], spherical, p),
    )
