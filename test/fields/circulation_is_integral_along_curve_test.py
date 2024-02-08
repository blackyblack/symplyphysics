from collections import namedtuple
from pytest import fixture
from sympy import sin, cos, pi
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
)
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.laws.fields import circulation_is_integral_along_curve as circulation_def


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    force_unit = Quantity(1 * units.newton)
    radius_unit = Quantity(1 * units.meter)
    # field is a field of gravitational forces, force is directed down by the Y coordinate
    # field is (0, -1 * G * m * M / y**2)
    # G * m * M = (force * length**2 / mass**2) * mass**2 = force * length**2
    field = VectorField(lambda point: [0, -1 * force_unit * radius_unit**2 / point.y**2], C)
    Args = namedtuple("Args", ["C", "radius_unit", "field"])
    return Args(C=C, radius_unit=radius_unit, field=field)


def test_basic_circulation(test_args):
    field = VectorField(lambda point: [point.y, 0, point.x + point.z], test_args.C)
    curve = [cos(circulation_def.parameter), sin(circulation_def.parameter)]
    result = circulation_def.calculate_circulation(field, curve, (0, pi / 2))
    assert_equal(result, -pi / 4)


def test_force_circulation(test_args):
    # trajectory is linear: y = x
    #HACK: gravitational force is undefined at 0 distance, use any non-zero value
    trajectory = [circulation_def.parameter, circulation_def.parameter]
    result = circulation_def.calculate_circulation(test_args.field, trajectory,
        (1 * test_args.radius_unit, 2 * test_args.radius_unit))
    assert_equal(result, -0.5 * units.newton * units.meter)


def test_force_circulation_horizontal(test_args):
    # trajectory is horizontal line: y = 5
    trajectory_horizontal = [circulation_def.parameter, 5 * test_args.radius_unit]
    result = circulation_def.calculate_circulation(test_args.field, trajectory_horizontal,
        (1 * test_args.radius_unit, 2 * test_args.radius_unit))
    assert_equal(result, 0)


def test_force_circulation_horizontal_up(test_args):
    # trajectory is vertical line: x = 5
    trajectory_vertical = [5 * test_args.radius_unit, circulation_def.parameter]
    result = circulation_def.calculate_circulation(test_args.field, trajectory_vertical,
        (1 * test_args.radius_unit, 2 * test_args.radius_unit))
    assert_equal(result, -0.5 * units.joule)


def test_force_circulation_horizontal_down(test_args):
    # trajectory is vertical line, but with down direction: x = 6
    trajectory_vertical = [6 * test_args.radius_unit, circulation_def.parameter]
    result = circulation_def.calculate_circulation(test_args.field, trajectory_vertical,
        (2 * test_args.radius_unit, 1 * test_args.radius_unit))
    assert_equal(result, 0.5 * units.joule)
