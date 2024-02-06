from collections import namedtuple
from pytest import fixture
from sympy import (sin, cos, pi, sqrt)
from symplyphysics import (
    assert_equal,
    SI,
    Quantity,
    units,
)
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.laws.fields import flux_is_integral_across_surface as flux_def


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    radius_unit = Quantity(1 * units.meter)
    Args = namedtuple("Args", ["C", "radius_unit"])
    return Args(C=C, radius_unit=radius_unit)


def test_basic_flux(test_args):
    field = VectorField(lambda point: [point.x, point.y, 0], test_args.C)
    # flux over hyperboloid, between planes z = -2 and z = 1
    surface = [
        sqrt(flux_def.parameter1**2 + 1) * cos(flux_def.parameter2),
        sqrt(flux_def.parameter1**2 + 1) * sin(flux_def.parameter2), -flux_def.parameter1
    ]
    result = flux_def.calculate_flux(field, surface, (-2, 1), (0, 2 * pi))
    assert_equal(result, 12 * pi)


def test_electric_intensity_flux(test_args):
    intensity = 100 * units.newton / units.coulomb
    field = VectorField([0, 2 * intensity, intensity], test_args.C)
    surface = [
        flux_def.parameter1 * cos(flux_def.parameter2),
        flux_def.parameter1 * sin(flux_def.parameter2), 0
    ]
    result = flux_def.calculate_flux(field, surface, (0, 0.02 * test_args.radius_unit), (0, 2 * pi))
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
        units.force * units.area / units.charge)
    assert_equal(result, (0.04 * pi) * units.newton * units.meter**2 / units.coulomb)
