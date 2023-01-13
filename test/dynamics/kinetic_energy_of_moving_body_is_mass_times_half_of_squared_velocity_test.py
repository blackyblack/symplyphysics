from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.dynamics import kinetic_energy_of_moving_body_is_mass_times_half_of_squared_velocity as kinetic_energy_law

# Description
## If we have a stone with 2kg of mass flying with 3m/s velocity, this stone should bear 9 Joules of kinetic energy.

@fixture
def test_args():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 2 * units.kilogram)

    V = units.Quantity('V')
    SI.set_quantity_dimension(V, units.velocity)
    SI.set_quantity_scale_factor(V, 3 * units.meter / units.second)

    Args = namedtuple('Args', ['m', 'V'])
    return Args(m = m, V = V)

def test_basic_energy(test_args):
    result = kinetic_energy_law.calculate_energy(test_args.m, test_args.V)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)

    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(2)
    assert result_energy == approx(9.0, 0.01)

def test_bad_mass(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 1 * units.meter)

    with raises(errors.UnitsError):
        kinetic_energy_law.calculate_energy(mb, test_args.V)

    with raises(TypeError):
        kinetic_energy_law.calculate_energy(100, test_args.V)

def test_bad_velocity(test_args):
    Vb = units.Quantity('Vb')
    SI.set_quantity_dimension(Vb, units.length)
    SI.set_quantity_scale_factor(Vb, 1 * units.meter)

    with raises(errors.UnitsError):
        kinetic_energy_law.calculate_energy(test_args.m, Vb)

    with raises(TypeError):
        kinetic_energy_law.calculate_energy(test_args.m, 100)