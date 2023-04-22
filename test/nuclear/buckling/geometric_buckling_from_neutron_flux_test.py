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
    sphere_radius = Quantity(10 * units.centimeter)
    neutron_flux = sin(pi / sphere_radius * spherical_coordinates.r) / spherical_coordinates.r
    Args = namedtuple("Args", ["f"])
    return Args(f=neutron_flux)


def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.f)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**-2)
    result_geometric_buckling = convert_to(result, units.centimeter**-2).subs(units.centimeter,
        1).evalf(2)
    assert result_geometric_buckling == approx(0.0986, 0.01)
