from collections import namedtuple
from pytest import approx, fixture
from math import pi
from sympy import sin, cos, sqrt
from symplyphysics import (
    units,
    convert_to,
    Quantity,
    SI,
    expr_to_quantity,
)
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.definitions import circulation_is_integral_along_curve as circulation_def


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    force_unit = Quantity(1 * units.newton)
    radius_unit = Quantity(1 * units.meter)
    # field is a field of gravitational forces, force is directed down by the Y coordinate
    # field is (0, -1 * G * m * M / y**2)
    # G * m * M = force * length**2 / mass**2 * mass**2 = force * length**2
    field = VectorField(0, lambda point: -1 * force_unit * radius_unit**2 / point.y**2, 0, C)
    Args = namedtuple("Args", ["C", "force_unit", "radius_unit", "field"])
    return Args(C=C, force_unit=force_unit, radius_unit=radius_unit, field=field)


def test_basic_circulation(test_args):
    field = VectorField(lambda point: point.y, 0, lambda point: point.x + point.z, test_args.C)
    curve = [cos(circulation_def.parameter), sin(circulation_def.parameter)]
    result_expr = circulation_def.calculate_circulation(field, curve, 0, pi / 2)
    assert result_expr.evalf(4) == approx(-pi / 4, 0.001)


def test_two_parameters_circulation(test_args):
    field = VectorField(lambda point: point.y, lambda point: -point.x, 0, test_args.C)
    # circle function is: x**2 + y**2 = 9
    # parametrize by circulation_def.parameter
    circle = [3 * cos(circulation_def.parameter), 3 * sin(circulation_def.parameter)]
    result_expr = circulation_def.calculate_circulation(field, circle, 0, 2 * pi)
    assert result_expr.evalf(4) == approx(-18 * pi, 0.001)
    # now try to define trajectory without parametrization
    # parametrized solution uses angle [0, 2*pi] that corresponds to the counter-clockwise direction
    # so we should integrate in the same direction: [r, -r] for upper part of the circle and [-r, r] for lower
    # y = sqrt(9 - x**2) for upper part of the circle
    # y = -sqrt(9 - x**2) for lower part of the circle
    circle_implicit_up = [circulation_def.parameter, sqrt(9 - circulation_def.parameter**2)]
    result_expr_up = circulation_def.calculate_circulation(field, circle_implicit_up, 3, -3)
    circle_implicit_down = [circulation_def.parameter, -sqrt(9 - circulation_def.parameter**2)]
    result_expr_down = circulation_def.calculate_circulation(field, circle_implicit_down, -3, 3)
    assert (result_expr_up + result_expr_down).evalf(4) == approx(-18 * pi, 0.001)


def test_orthogonal_movement_circulation(test_args):
    field = VectorField(lambda point: point.y, lambda point: -point.x, 1, test_args.C)
    # trajectory is upwards helix
    helix = [
        cos(circulation_def.parameter),
        sin(circulation_def.parameter), circulation_def.parameter
    ]
    result_expr = circulation_def.calculate_circulation(field, helix, 0, 2 * pi)
    assert result_expr == 0
    # trajectory is upwards straight line
    trajectory_vertical = [1, 0, circulation_def.parameter]
    result_expr = circulation_def.calculate_circulation(field, trajectory_vertical, 0, 2 * pi)
    assert result_expr.evalf(4) == approx(2 * pi, 0.001)


def test_force_circulation(test_args):
    # trajectory is linear: y = x
    #HACK: gravitational force is undefined at 0 distance, use any non-zero value
    trajectory = [circulation_def.parameter, circulation_def.parameter]
    result_expr = circulation_def.calculate_circulation(test_args.field, trajectory,
        1 * test_args.radius_unit, 2 * test_args.radius_unit)
    result = expr_to_quantity(result_expr)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).evalf(2)
    assert result_work == approx(-0.5, 0.01)


def test_force_circulation_horizontal(test_args):
    # trajectory is horizontal line: y = 5
    trajectory_horizontal = [circulation_def.parameter, 5 * test_args.radius_unit]
    result_expr = circulation_def.calculate_circulation(test_args.field, trajectory_horizontal,
        1 * test_args.radius_unit, 2 * test_args.radius_unit)
    assert result_expr == 0


def test_force_circulation_horizontal_up(test_args):
    # trajectory is vertical line: x = 5
    trajectory_vertical = [5 * test_args.radius_unit, circulation_def.parameter]
    result_expr = circulation_def.calculate_circulation(test_args.field, trajectory_vertical,
        1 * test_args.radius_unit, 2 * test_args.radius_unit)
    result = expr_to_quantity(result_expr)
    result_work = convert_to(result, units.joule).evalf(2)
    assert result_work == approx(-0.5, 0.01)


def test_force_circulation_horizontal_down(test_args):
    # trajectory is vertical line, but with down direction: x = 6
    trajectory_vertical = [6 * test_args.radius_unit, circulation_def.parameter]
    result_expr = circulation_def.calculate_circulation(test_args.field, trajectory_vertical,
        2 * test_args.radius_unit, 1 * test_args.radius_unit)
    result = expr_to_quantity(result_expr)
    result_work = convert_to(result, units.joule).evalf(2)
    assert result_work == approx(0.5, 0.01)
