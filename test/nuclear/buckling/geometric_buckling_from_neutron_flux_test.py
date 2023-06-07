from collections import namedtuple
from pytest import approx, fixture
from sympy import sin, pi
from sympy.vector import CoordSys3D
from symplyphysics import (
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux as buckling


@fixture
def test_args():
    # spherical reactor with radius = 10 centimeter
    spherical_coordinates = CoordSys3D("spherical_coordinates", transformation="spherical")
    unit_length = Quantity(1 * units.meter)
    # neutron flux function has radius in denominator, hence the coefficient should be multiplied to units.length
    # to conform to neutron flux dimension
    neutron_flux_times_radius_unit = Quantity(1 / units.meter / units.second)
    sphere_radius = Quantity(10 * units.centimeter)
    neutron_flux = neutron_flux_times_radius_unit * sin(pi / sphere_radius *
        spherical_coordinates.r * unit_length) / (spherical_coordinates.r * unit_length)
    Args = namedtuple("Args", ["f"])
    return Args(f=neutron_flux)


def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.f)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.length**2)
    result_geometric_buckling = convert_to(result,
        1 / units.centimeter**2).subs(units.centimeter, 1).evalf(4)
    assert result_geometric_buckling == approx(0.0986, 0.01)
