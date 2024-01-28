from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (units, SI, convert_to, Quantity, errors)

from symplyphysics.laws.electricity import force_from_charge_velocity_magnetic_induction as force_law

# Description
## The charge value is 0.1 coulomb. The charge rate is 3 meter per second. The magnetic induction is 2 tesla.
## The angle between induction and velocity is 45 degree (pi / 4 radian). Then the Lorentz force is 0.4 newton.
## https://physics.icalculator.com/lorentz-force-calculator.html


@fixture(name="test_args")
def test_args_fixture():
    charge = Quantity(0.1 * units.coulomb)
    velocity = Quantity(3 * (units.meter / units.second))
    induction = Quantity(2 * units.tesla)
    angle = pi / 4

    Args = namedtuple("Args", ["charge", "velocity", "induction", "angle"])
    return Args(charge=charge, velocity=velocity, angle=angle, induction=induction)


def test_basic_force(test_args):
    result = force_law.calculate_force(test_args.charge, test_args.velocity, test_args.angle,
        test_args.induction)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result = convert_to(result, units.newton).evalf(5)
    assert result == approx(0.4, rel=0.1)


def test_bad_charge(test_args):
    charge = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        force_law.calculate_force(charge, test_args.velocity, test_args.angle, test_args.induction)
    with raises(TypeError):
        force_law.calculate_force(100, test_args.velocity, test_args.angle, test_args.induction)


def test_bad_velocity(test_args):
    velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.charge, velocity, test_args.angle, test_args.induction)
    with raises(TypeError):
        force_law.calculate_force(test_args.charge, 100, test_args.angle, test_args.induction)


def test_bad_angle(test_args):
    angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.charge, test_args.velocity, angle, test_args.induction)
    with raises(AttributeError):
        force_law.calculate_force(test_args.charge, test_args.velocity, True, test_args.induction)


def test_bad_induction(test_args):
    induction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.charge, test_args.velocity, test_args.angle, induction)
    with raises(TypeError):
        force_law.calculate_force(test_args.charge, test_args.velocity, test_args.angle, 100)
