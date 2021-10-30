from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.chemistry import atomic_weight_from_mass_mole_count

@fixture
def test_args():
    mass = units.Quantity('mass')
    SI.set_quantity_dimension(mass, units.mass)
    # molar mass of water is 18.0153 gram / mole
    SI.set_quantity_scale_factor(mass, 18.0153 * units.gram)
    mole_count = units.Quantity('mole_count')
    SI.set_quantity_dimension(mole_count, units.amount_of_substance)
    SI.set_quantity_scale_factor(mole_count, 1 * units.mole)

    Args = namedtuple('Args', ['m', 'N'])
    return Args(m = mass, N = mole_count)

def test_basic_atomic_weight(test_args):
    result = atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, test_args.N)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass / units.amount_of_substance)

    result_mass = convert_to(result, units.gram / units.mole).subs(units.gram, 1).subs(units.mole, 1).evalf(4)
    assert result_mass == approx(18.0153, 0.01)

def test_bad_mass(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 3 * units.meter)

    with raises(errors.UnitsError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(mb, test_args.N)

    with raises(TypeError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(100, test_args.N)

def test_bad_mole_count(test_args):
    Nb = units.Quantity('Nb')
    SI.set_quantity_dimension(Nb, units.length)
    SI.set_quantity_scale_factor(Nb, 3 * units.meter)

    with raises(errors.UnitsError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, Nb)

    with raises(TypeError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, 100)
