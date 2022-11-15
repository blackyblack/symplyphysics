from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.kinematics import kinetic_energy_of_body_from_mass_and_velocity as kinetic_energy

@fixture
def test_args():
    body_mass = units.Quantity('body_mass')
    SI.set_quantity_dimension(body_mass, units.mass)
    SI.set_quantity_scale_factor(body_mass, 0.5 * units.kilogram)
    body_velocity = units.Quantity('body_velocity')
    SI.set_quantity_dimension(body_velocity, units.velocity)
    SI.set_quantity_scale_factor(body_velocity, 0.5 * units.meter / units.second)
    Args = namedtuple('Args', ['body_mass', 'body_velocity'])
    return Args(body_mass=body_mass, body_velocity=body_velocity)

def test_basic_kinematic_energy(test_args):
    result = kinetic_energy.calculate_kinetic_energy_body(test_args.body_mass, test_args.body_velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(3)
    assert result_energy == approx(0.0625, 0.005)

def test_bad_body_mass(test_args):
    bbm = units.Quantity('bbm')
    SI.set_quantity_dimension(bbm, units.velocity)
    SI.set_quantity_scale_factor(bbm, 1 * units.meter)
    with raises(errors.UnitsError):
        kinetic_energy.calculate_kinetic_energy_body(bbm, test_args.body_velocity)
    with raises(TypeError):
        kinetic_energy.calculate_kinetic_energy_body(100, test_args.body_velocity)

def test_bad_body_velocity(test_args):
    bbv = units.Quantity('bbv')
    SI.set_quantity_dimension(bbv, units.mass)
    SI.set_quantity_scale_factor(bbv, 1 * units.meter)
    with raises(errors.UnitsError):
        kinetic_energy.calculate_kinetic_energy_body(test_args.body_mass, bbv)
    with raises(TypeError):
        kinetic_energy.calculate_kinetic_energy_body(test_args.body_mass, 100)
