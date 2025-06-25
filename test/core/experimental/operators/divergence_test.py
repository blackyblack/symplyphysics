from pytest import raises

from sympy import sin, cos, tan, Function as SymFunction

from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.symbols import BasicSymbol, Symbol

from symplyphysics.core.experimental.vectors import VectorSymbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    CylindricalCoordinateSystem, SphericalCoordinateSystem)
from symplyphysics.core.experimental.coordinate_systems.vector import CoordinateVector
from symplyphysics.core.experimental.coordinate_systems.scalar import CoordinateScalar
from symplyphysics.core.experimental.operators import VectorDivergence
from symplyphysics.core.experimental.solvers import vector_equals


def test_non_coordinate_vector() -> None:
    assert vector_equals(VectorDivergence(0), 0)

    with raises(ValueError):
        assert vector_equals(VectorDivergence(1), 0)

    x = Symbol("x", real=True)
    with raises(ValueError):
        _ = VectorDivergence(x)

    f = SymFunction("f")
    with raises(ValueError):
        _ = VectorDivergence(f(x))

    v = VectorSymbol("v")
    assert VectorDivergence(v) == 0

    w = VectorSymbol("w")
    assert VectorDivergence(v + w) == 0

    assert VectorDivergence(VectorCross(v, w)) == 0


def test_non_evaluate() -> None:
    assert vector_equals(VectorDivergence(0, evaluate=False), 0)

    v = VectorSymbol("v")
    div = VectorDivergence(v, evaluate=False)
    assert isinstance(div, VectorDivergence)
    assert div.args == (v,)
    assert expr_equals(div, 0)  # evaluate

    w = VectorSymbol("w")
    # v_plus_w = v + w
    div = VectorDivergence(v + w, evaluate=False)
    assert isinstance(div, VectorDivergence)
    assert div.args == (v + w,)
    assert expr_equals(div, 0)  # evaluate

    div = VectorDivergence(VectorCross(v, w), evaluate=False)
    assert isinstance(div, VectorDivergence)
    assert div.args == (VectorCross(v, w),)
    assert expr_equals(div, 0)  # evaluate

    assert VectorDivergence(VectorCross(v, w)) == 0

    cartesian = CartesianCoordinateSystem()
    x, y, z = cartesian.base_scalars
    vector = CoordinateVector([x, y, z], cartesian)
    div = VectorDivergence(vector, evaluate=False)
    assert isinstance(div, VectorDivergence)
    assert div.args == (vector,)


def test_cartesian_system() -> None:
    cartesian = CartesianCoordinateSystem()

    x, y, z = cartesian.base_scalars

    vector: CoordinateVector = CoordinateVector([x, y, -z], cartesian)
    div = VectorDivergence(vector)
    assert isinstance(div, CoordinateScalar)
    assert expr_equals(div.scalar, 1)
    assert div.system == cartesian
    assert div.point == vector.point

    vector = CoordinateVector([y, z, x], cartesian)
    assert VectorDivergence(vector) == 0

    vector = CoordinateVector([x * y, y * z, z * x], cartesian)
    div = VectorDivergence(vector)
    assert isinstance(div, CoordinateScalar)
    assert expr_equals(div.scalar, x + y + z)
    assert div.system == cartesian
    assert div.point == vector.point


def test_cylindrical_system() -> None:
    cylindrical = CylindricalCoordinateSystem()

    rho, phi, z = cylindrical.base_scalars
    p = BasicSymbol("P")

    vector: CoordinateVector = CoordinateVector([rho, z * sin(phi), z * cos(phi)], cylindrical, p)
    div = VectorDivergence(vector)
    assert isinstance(div, CoordinateScalar)
    assert expr_equals(div.scalar, 2 + cos(phi) + z * cos(phi) / rho)
    assert div.system == cylindrical
    assert div.point == vector.point

    vector = CoordinateVector([-1 * rho * cos(phi) / 2, rho * sin(phi), rho], cylindrical, p)
    assert VectorDivergence(vector) == 0


def test_spherical_system() -> None:
    spherical = SphericalCoordinateSystem()

    r, theta, phi = spherical.base_scalars
    p = BasicSymbol("P")

    vector: CoordinateVector = CoordinateVector([r, r * cos(theta), r * sin(phi)], spherical, p)
    div = VectorDivergence(vector)
    assert isinstance(div, CoordinateScalar)
    assert expr_equals(div.scalar, 3 + (cos(phi) + cos(2 * theta)) / sin(theta))
    assert div.system == spherical
    assert div.point == vector.point

    vector = CoordinateVector([r / 3, r / tan(theta), r * sin(theta)], spherical, p)
    div = VectorDivergence(vector)
    assert isinstance(div, CoordinateScalar)
    assert expr_equals(div.scalar, 0)
    assert div.system == spherical
    assert div.point == vector.point
