from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import radius_of_a_charged_particle_in_a_magnetic_field as radius_law

# Description
## The mass of the particle is 0.1 kilogram, the charge is 0.2 coulomb, the speed is 5 meter per second.
## With a magnetic induction equal to 3 tesla, the radius of motion of the particle is 0.8333 meter.
## https://physics.icalculator.com/radius-of-trajectory-and-period-of-a-charge-moving-inside-a-uniform-magnetic-field-calculator.html

Args = namedtuple("Args", ["mass", "velocity", "charge", "induction"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass = Quantity(0.1 * units.kilogram)
    velocity = Quantity(5 * (units.meter / units.second))
    charge = Quantity(0.2 * units.coulomb)
    induction = Quantity(3 * units.tesla)
    return Args(mass=mass, velocity=velocity, induction=induction, charge=charge)


def test_basic_radius(test_args: Args) -> None:
    result = radius_law.calculate_radius(test_args.mass, test_args.velocity, test_args.induction,
        test_args.charge)
    assert_equal(result, 0.8333 * units.meter)


def test_bad_mass(test_args: Args) -> None:
    mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(mass, test_args.velocity, test_args.induction, test_args.charge)
    with raises(TypeError):
        radius_law.calculate_radius(100, test_args.velocity, test_args.induction, test_args.charge)


def test_bad_velocity(test_args: Args) -> None:
    velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.mass, velocity, test_args.induction, test_args.charge)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.mass, 100, test_args.induction, test_args.charge)


def test_bad_induction(test_args: Args) -> None:
    induction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.mass, test_args.velocity, induction, test_args.charge)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.mass, test_args.velocity, 100, test_args.charge)


def test_bad_charge(test_args: Args) -> None:
    charge = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.mass, test_args.velocity, test_args.induction, charge)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.mass, test_args.velocity, test_args.induction, 100)
