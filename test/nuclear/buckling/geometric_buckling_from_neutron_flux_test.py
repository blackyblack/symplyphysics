from collections import namedtuple
from pytest import fixture
from sympy import sin, pi
from sympy.vector import CoordSys3D
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux as buckling


@fixture(name="test_args")
def test_args_fixture():
    # spherical reactor with radius = 10 centimeter
    spherical_coordinates = CoordSys3D("spherical_coordinates", transformation="spherical")
    # Makes linter happy
    r = getattr(spherical_coordinates, "r")
    unit_length = Quantity(1 * units.meter)
    # neutron flux function has radius in denominator, hence the coefficient should be multiplied to units.length
    # to conform to neutron flux dimension
    neutron_flux_times_radius_unit = Quantity(1 / units.meter / units.second)
    sphere_radius = Quantity(10 * units.centimeter)
    neutron_flux = neutron_flux_times_radius_unit * sin(
        pi / sphere_radius * r * unit_length) / (r * unit_length)
    Args = namedtuple("Args", ["f"])
    return Args(f=neutron_flux)


def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.f)
    assert_equal(result, 0.0986 / units.centimeter**2)
