from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)

from symplyphysics.laws.electricity import force_from_charge_velocity_magnetic_induction as force_law

# Description
## The charge value is 0.1 coulomb. The charge rate is 3 meter per second. The magnetic induction is 2 tesla.
## The angle between induction and velocity is 45 degree (pi / 4 radian). Then the Lorentz force is 0.424 newton.
## https://physics.icalculator.com/lorentz-force-calculator.html

Args = namedtuple("Args", ["charge", "velocity", "induction", "angle"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    charge = Quantity(0.1 * units.coulomb)
    velocity = Quantity(3 * (units.meter / units.second))
    induction = Quantity(2 * units.tesla)
    angle = pi / 4
    return Args(charge=charge, velocity=velocity, angle=angle, induction=induction)


def test_basic_force(test_args: Args) -> None:
    result = force_law.calculate_force(test_args.charge, test_args.velocity, test_args.angle,
        test_args.induction)
    assert_equal(result, 0.424 * units.newton)


def test_bad_charge(test_args: Args) -> None:
    charge = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        force_law.calculate_force(charge, test_args.velocity, test_args.angle, test_args.induction)
    with raises(TypeError):
        force_law.calculate_force(100, test_args.velocity, test_args.angle, test_args.induction)


def test_bad_velocity(test_args: Args) -> None:
    velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.charge, velocity, test_args.angle, test_args.induction)
    with raises(TypeError):
        force_law.calculate_force(test_args.charge, 100, test_args.angle, test_args.induction)


def test_bad_angle(test_args: Args) -> None:
    angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.charge, test_args.velocity, angle, test_args.induction)
    with raises(AttributeError):
        force_law.calculate_force(test_args.charge, test_args.velocity, True, test_args.induction)


def test_bad_induction(test_args: Args) -> None:
    induction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.charge, test_args.velocity, test_args.angle, induction)
    with raises(TypeError):
        force_law.calculate_force(test_args.charge, test_args.velocity, test_args.angle, 100)
