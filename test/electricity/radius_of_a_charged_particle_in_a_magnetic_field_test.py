from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.electricity import radius_of_a_charged_particle_in_a_magnetic_field as radius_law

# Description
## The mass of the particle is 0.1 kilogram, the charge is 0.2 coulomb, the speed is 5 meter per second.
## With a magnetic induction equal to 3 tesla, the radius of motion of the particle is 0.83 meter.
## https://physics.icalculator.com/radius-of-trajectory-and-period-of-a-charge-moving-inside-a-uniform-magnetic-field-calculator.html


@fixture(name="test_args")
def test_args_fixture():
    mass = Quantity(0.1 * units.kilogram)
    velocity = Quantity(5 * (units.meter / units.second))
    charge = Quantity(0.2 * units.coulomb)
    induction = Quantity(3 * units.tesla)

    Args = namedtuple("Args", ["mass", "velocity", "charge", "induction"])
    return Args(mass=mass, velocity=velocity, induction=induction, charge=charge)


def test_basic_radius(test_args):
    result = radius_law.calculate_radius(test_args.mass, test_args.velocity, test_args.induction,
        test_args.charge)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result = convert_to(result, units.meter).evalf(5)
    assert result == approx(0.83, rel=0.1)


def test_bad_mass(test_args):
    mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(mass, test_args.velocity, test_args.induction, test_args.charge)
    with raises(TypeError):
        radius_law.calculate_radius(100, test_args.velocity, test_args.induction, test_args.charge)


def test_bad_velocity(test_args):
    velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.mass, velocity, test_args.induction, test_args.charge)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.mass, 100, test_args.induction, test_args.charge)


def test_bad_induction(test_args):
    induction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.mass, test_args.velocity, induction, test_args.charge)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.mass, test_args.velocity, 100, test_args.charge)


def test_bad_charge(test_args):
    charge = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.mass, test_args.velocity, test_args.induction, charge)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.mass, test_args.velocity, test_args.induction, 100)
