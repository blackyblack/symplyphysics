from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import density_from_mass_volume

@fixture
def test_args():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 1 * units.kilogram)
    V = units.Quantity('V')
    SI.set_quantity_dimension(V, units.volume)
    SI.set_quantity_scale_factor(V, 3 * units.meter**3)

    Args = namedtuple('Args', ['m', 'V'])
    return Args(m=m, V=V)


def test_basic_density(test_args):
    result = density_from_mass_volume.calculate_density(test_args.m, test_args.V)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass / units.volume)

    result_density = convert_to(result, density_from_mass_volume.definition_dimension_SI).subs({
        units.kilogram: 1, units.meter: 1}).evalf(2)
    assert result_density == approx(0.3333, 0.01)

def test_bad_mass(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 1 * units.meter)
    with raises(errors.UnitsError):
        density_from_mass_volume.calculate_density(mb, test_args.V)

    with raises(TypeError):
        density_from_mass_volume.calculate_density(100, test_args.V)

def test_bad_volume(test_args):
    Vb = units.Quantity('Vb')
    SI.set_quantity_dimension(Vb, units.length)
    SI.set_quantity_scale_factor(Vb, 1 * units.meter)
    with raises(errors.UnitsError):
        density_from_mass_volume.calculate_density(test_args.m, Vb)

    with raises(TypeError):
        density_from_mass_volume.calculate_density(test_args.m, 100)
