from collections import namedtuple
from pytest import fixture
from sympy import Expr, sin, cos, sqrt, pi
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    CoordinateSystem,
    SI,
    convert_to,
)
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.points.cartesian_point import CartesianPoint
from symplyphysics.laws.fields import circulation_is_integral_of_curl_over_surface as circulation_def


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    force_unit = Quantity(1 * units.newton)
    radius_unit = Quantity(1 * units.meter)
    Args = namedtuple("Args", ["C", "force_unit", "radius_unit"])
    return Args(C=C, force_unit=force_unit, radius_unit=radius_unit)


def test_basic_circulation(test_args):
    field = VectorField(lambda point: [point.y, 0, point.x + point.z], test_args.C)
    surface = [
        circulation_def.parameter1 * cos(circulation_def.parameter2),
        circulation_def.parameter1 * sin(circulation_def.parameter2)
    ]
    result = circulation_def.calculate_circulation(field, surface, (0, 1), (0, pi / 2))
    assert_equal(result, -pi / 4)


def _distance(point: CartesianPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)


def test_force_field_circulation(test_args):
    # we use lorentz force in magnetic field as reference
    # B = mass / (current * time**2) = mass / (charge * time)
    # Lorentz force is: F = q * v * B = charge * (length / time) * B
    field = VectorField(
        lambda point: [
        -point.y * test_args.force_unit / _distance(point), point.x * test_args.force_unit /
        _distance(point), 0
        ], test_args.C)
    surface = [
        circulation_def.parameter1 * cos(circulation_def.parameter2),
        circulation_def.parameter1 * sin(circulation_def.parameter2)
    ]
    result = circulation_def.calculate_circulation(field, surface,
        (1 * test_args.radius_unit, 2 * test_args.radius_unit), (0, pi / 2))
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).evalf(2)
    assert result_work > 0
