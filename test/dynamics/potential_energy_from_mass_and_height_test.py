from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy

# How much potential energy does the body with a mass of 9 grams have at
# an altitude of 500 meters?

@fixture
def test_args():
    m = Quantity(units.mass, 9 * units.gram)
    h = Quantity(units.length, 500 * units.meter)
    Args = namedtuple("Args", ["m", "h"])
    return Args(m=m, h=h)

def test_basic_energy(test_args):
    result = potential_energy.calculate_potential_energy(test_args.m, test_args.h)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(3)
    assert result_energy == approx(44.19, 0.01)

def test_bad_body_mass(test_args):
    bm = Quantity(units.length)
    with raises(errors.UnitsError):
        potential_energy.calculate_potential_energy(bm, test_args.h)
    with raises(TypeError):
        potential_energy.calculate_potential_energy(100, test_args.h)

def test_bad_height(test_args):
    bh = Quantity(units.mass)
    with raises(errors.UnitsError):
        potential_energy.calculate_potential_energy(test_args.m, bh)
    with raises(TypeError):
        potential_energy.calculate_potential_energy(test_args.m, 100)
