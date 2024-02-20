from sympy import Expr, symbols
from sympy.vector import CoordSys3D
from symplyphysics.core.points.cartesian_point import CartesianPoint
from symplyphysics.core.points.cylinder_point import CylinderPoint
from symplyphysics.core.points.sphere_point import SpherePoint
from symplyphysics.core.points.point import Point


def _assert_point(point_: Point, expected_: list[Expr | float]) -> None:
    for idx, c in enumerate(point_.coordinates):
        assert c == expected_[idx]


def test_basic_point() -> None:
    point = Point(1, 2, 3)
    _assert_point(point, [1, 2, 3])
    assert point.coordinate(0) == 1
    assert point.coordinate(1) == 2
    assert point.coordinate(2) == 3


def test_basic_cartesian_point() -> None:
    point = CartesianPoint(1, 2, 3)
    _assert_point(point, [1, 2, 3])
    assert point.x == 1
    assert point.y == 2
    assert point.z == 3


def test_basic_sphere_point() -> None:
    point = SpherePoint(1, 2, 3)
    _assert_point(point, [1, 2, 3])
    assert point.r == 1
    assert point.theta == 2
    assert point.phi == 3
    assert point.radius == 1
    assert point.azimuthal_angle == 2
    assert point.polar_angle == 3


def test_basic_cylinder_point() -> None:
    point = CylinderPoint(1, 2, 3)
    _assert_point(point, [1, 2, 3])
    assert point.r == 1
    assert point.theta == 2
    assert point.z == 3
    assert point.radius == 1
    assert point.azimuthal_angle == 2
    assert point.height == 3


def test_empty_point() -> None:
    point = Point()
    _assert_point(point, [0, 0, 0])


def test_4d_point() -> None:
    point = Point(1, 2, 3, 4)
    _assert_point(point, [1, 2, 3, 4])


# This kind of FieldPoint represents a surface, so it can contain any point in XY-plane
# and it's Z-coordinate = 3
def test_coord_system_point() -> None:
    C = CoordSys3D("C")
    # Make linter happy
    x = getattr(C, "x")
    y = getattr(C, "y")
    point = Point(x, y, 3)
    _assert_point(point, [x, y, 3])


# This kind of FieldPoint represents a straight line in 2D space. It's function is Y = X,
# as X = param and Y = param. It's Z-coordinate = 0.
def test_parametrized_point() -> None:
    param = symbols("param")
    point = Point(param, param)
    _assert_point(point, [param, param, 0])
