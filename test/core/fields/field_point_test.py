from sympy import Expr, symbols
from sympy.vector import CoordSys3D
from symplyphysics.core.fields.field_point import FieldPoint


def _assert_point(point_: FieldPoint, expected_: list[Expr | float]):
    for idx, c in enumerate(point_.coordinates):
        assert c == expected_[idx]


def test_basic_point():
    point = FieldPoint(1, 2, 3)
    _assert_point(point, [1, 2, 3])
    assert point.x == 1
    assert point.y == 2
    assert point.z == 3


def test_empty_point():
    point = FieldPoint()
    _assert_point(point, [0, 0, 0])


def test_skip_coordinate_point():
    point = FieldPoint(1, None, 3)
    _assert_point(point, [1, 0, 3])


def test_4d_point():
    point = FieldPoint(1, 2, 3)
    point.set_coordinate(3, 4)
    _assert_point(point, [1, 2, 3, 4])


# This kind of FieldPoint represents a surface, so it can contain any point in XY-plane
# and it's Z-coordinate = 3
def test_coord_system_point():
    C = CoordSys3D("C")
    # Make linter happy
    x = getattr(C, "x")
    y = getattr(C, "y")
    point = FieldPoint(x, y, 3)
    _assert_point(point, [x, y, 3])


# This kind of FieldPoint represents a straight line in 2D space. It's function is Y = X,
# as X = param and Y = param. It's Z-coordinate = 0.
def test_parametrized_point():
    param = symbols("param")
    point = FieldPoint(param, param)
    _assert_point(point, [param, param, 0])
