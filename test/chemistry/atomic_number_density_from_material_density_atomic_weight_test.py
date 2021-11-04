from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.chemistry import atomic_number_density_from_material_density_atomic_weight as atomic_number_density

@fixture
def test_args():
    p = units.Quantity('p')
    SI.set_quantity_dimension(p, units.mass / units.volume)
    # boron carbide density is 2.52 g/cm^3
    SI.set_quantity_scale_factor(p, 2.52 * units.gram / units.centimeter**3)
    M = units.Quantity('M')
    SI.set_quantity_dimension(M, units.mass / units.amount_of_substance)
    # boron carbide atomic mass
    SI.set_quantity_scale_factor(M, 55.2 * units.gram / units.mole)

    Args = namedtuple('Args', ['p', 'M'])
    return Args(p = p, M = M)

def test_basic_number_density(test_args):
    result = atomic_number_density.calculate_atomic_number_density(test_args.p, test_args.M)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.volume)

    result_number_density = convert_to(result, 1 / units.centimeter**3).subs(units.centimeter, 1).evalf(4)
    assert result_number_density == approx(2.75E+22, 0.01)

def test_bad_density(test_args):
    pb = units.Quantity('pb')
    SI.set_quantity_dimension(pb, units.length)
    SI.set_quantity_scale_factor(pb, 1 * units.meter)

    with raises(errors.UnitsError):
        atomic_number_density.calculate_atomic_number_density(pb, test_args.M)

    with raises(TypeError):
        atomic_number_density.calculate_atomic_number_density(100, test_args.M)

def test_bad_atomic_mass(test_args):
    Mb = units.Quantity('Mb')
    SI.set_quantity_dimension(Mb, units.length)
    SI.set_quantity_scale_factor(Mb, 3 * units.meter)

    with raises(errors.UnitsError):
        atomic_number_density.calculate_atomic_number_density(test_args.p, Mb)

    with raises(TypeError):
        atomic_number_density.calculate_atomic_number_density(test_args.p, 100)
