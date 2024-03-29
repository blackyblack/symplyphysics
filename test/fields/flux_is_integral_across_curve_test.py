from collections import namedtuple
from pytest import fixture
from sympy import sin, cos, pi
from symplyphysics import (
    assert_equal,
    Quantity,
    units,
    SI,
)
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.laws.fields import flux_is_integral_across_curve as flux_def

Args = namedtuple("Args", ["C", "force_unit", "radius_unit"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = CoordinateSystem()
    force_unit = Quantity(1 * units.newton)
    radius_unit = Quantity(1 * units.meter)
    return Args(C=C, force_unit=force_unit, radius_unit=radius_unit)


def test_basic_flux(test_args: Args) -> None:
    field = VectorField(lambda point: [point.x, point.y], test_args.C)
    # flux over circle of radius = 2
    curve = [2 * cos(flux_def.parameter), 2 * sin(flux_def.parameter)]
    result = flux_def.calculate_flux(field, curve, (0, 2 * pi))
    assert_equal(result, 8 * pi)


def test_gravitational_field_flux(test_args: Args) -> None:
    field = VectorField(
        lambda point: [0, -1 * test_args.force_unit * test_args.radius_unit**2 / point.y**2],
        test_args.C)
    trajectory = [flux_def.parameter, flux_def.parameter]
    result = flux_def.calculate_flux(field, trajectory,
        (1 * test_args.radius_unit, 2 * test_args.radius_unit))
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
