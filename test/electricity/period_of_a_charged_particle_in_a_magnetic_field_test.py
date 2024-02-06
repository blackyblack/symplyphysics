from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import period_of_a_charged_particle_in_a_magnetic_field as period_law

# Description
## The mass of the particle is 0.1 kilogram, the charge is 0.2 coulomb.
## With a magnetic induction equal to 3 tesla, the period of motion of the particle is 1.047 second.
## https://physics.icalculator.com/radius-of-trajectory-and-period-of-a-charge-moving-inside-a-uniform-magnetic-field-calculator.html


@fixture(name="test_args")
def test_args_fixture():
    mass = Quantity(0.1 * units.kilogram)
    charge = Quantity(0.2 * units.coulomb)
    induction = Quantity(3 * units.tesla)

    Args = namedtuple("Args", ["mass", "charge", "induction"])
    return Args(mass=mass, charge=charge, induction=induction)


def test_basic_period(test_args):
    result = period_law.calculate_period(test_args.mass, test_args.charge, test_args.induction)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.time)
    result_period = convert_to(result, units.second).evalf(5)
    assert_approx(result_period, 1.047)


def test_bad_mass(test_args):
    mass = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        period_law.calculate_period(mass, test_args.charge, test_args.induction)
    with raises(TypeError):
        period_law.calculate_period(100, test_args.charge, test_args.induction)


def test_bad_charge(test_args):
    charge = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        period_law.calculate_period(test_args.mass, charge, test_args.induction)
    with raises(TypeError):
        period_law.calculate_period(test_args.mass, 100, test_args.induction)


def test_bad_induction(test_args):
    induction = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        period_law.calculate_period(test_args.mass, test_args.charge, induction)
    with raises(TypeError):
        period_law.calculate_period(test_args.mass, test_args.charge, 100)
