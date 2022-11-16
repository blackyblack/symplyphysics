from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.kinematic import kinetic_energy_from_mass_and_velocity as kinetic_energy

@fixture
def test_args():
    mass = units.Quantity('mass')
    SI.set_quantity_dimension(mass, units.mass)
    SI.set_quantity_scale_factor(mass, 0.5 * units.kilogram)
    velocity = units.Quantity('velocity')
    SI.set_quantity_dimension(velocity, units.velocity)
    SI.set_quantity_scale_factor(velocity, 0.5 * units.meter / units.second)
    Args = namedtuple('Args', ['mass', 'velocity'])
    return Args(mass=mass, velocity=velocity)

def test_basic_kinetic_energy(test_args):
    result = kinetic_energy.calculate_kinetic_energy(test_args.mass, test_args.velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(3)
    assert result_energy == approx(0.0625, 0.005)

def test_bad_body_mass(test_args):
    bbm = units.Quantity('bbm')
    SI.set_quantity_dimension(bbm, units.velocity)
    SI.set_quantity_scale_factor(bbm, 1 * units.meter)
    with raises(errors.UnitsError):
        kinetic_energy.calculate_kinetic_energy(bbm, test_args.velocity)
    with raises(TypeError):
        kinetic_energy.calculate_kinetic_energy(100, test_args.velocity)

def test_bad_body_velocity(test_args):
    bbv = units.Quantity('bbv')
    SI.set_quantity_dimension(bbv, units.mass)
    SI.set_quantity_scale_factor(bbv, 1 * units.meter)
    with raises(errors.UnitsError):
        kinetic_energy.calculate_kinetic_energy(test_args.mass, bbv)
    with raises(TypeError):
        kinetic_energy.calculate_kinetic_energy(test_args.mass, 100)
