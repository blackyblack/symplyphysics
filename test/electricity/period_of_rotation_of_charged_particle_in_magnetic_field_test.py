from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import period_of_rotation_of_charged_particle_in_magnetic_field as period_law

# Description
## The mass of the particle is 0.1 kilogram, the charge is 0.2 coulomb.
## With a magnetic induction equal to 3 tesla, the period of motion of the particle is 1.047 second.
## https://physics.icalculator.com/radius-of-trajectory-and-period-of-a-charge-moving-inside-a-uniform-magnetic-field-calculator.html

Args = namedtuple("Args", ["mass", "charge", "induction"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass = Quantity(0.1 * units.kilogram)
    charge = Quantity(0.2 * units.coulomb)
    induction = Quantity(3 * units.tesla)
    return Args(mass=mass, charge=charge, induction=induction)


def test_basic_period(test_args: Args) -> None:
    result = period_law.calculate_period(test_args.mass, test_args.charge, test_args.induction)
    assert_equal(result, 1.047 * units.second)


def test_bad_mass(test_args: Args) -> None:
    mass = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        period_law.calculate_period(mass, test_args.charge, test_args.induction)
    with raises(TypeError):
        period_law.calculate_period(100, test_args.charge, test_args.induction)


def test_bad_charge(test_args: Args) -> None:
    charge = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        period_law.calculate_period(test_args.mass, charge, test_args.induction)
    with raises(TypeError):
        period_law.calculate_period(test_args.mass, 100, test_args.induction)


def test_bad_induction(test_args: Args) -> None:
    induction = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        period_law.calculate_period(test_args.mass, test_args.charge, induction)
    with raises(TypeError):
        period_law.calculate_period(test_args.mass, test_args.charge, 100)
