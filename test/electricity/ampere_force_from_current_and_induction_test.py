from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import ampere_force_from_current_and_induction

# Description
## With a current of 3.5 amperes, a conductor length of 2 meters and a magnetic induction of 0.6061 tesla,
## the force will be equal to 3 newton. The angle between the magnetic induction and the current direction
## is 45 degrees (pi / 4 radians).
## https://physics.icalculator.com/calculating-magnetic-field-using-the-amperes-law.html

Args = namedtuple("Args", ["current", "length", "induction", "angle"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    current = Quantity(3.5 * units.ampere)
    length = Quantity(2 * units.meter)
    induction = Quantity(0.6061 * units.tesla)
    angle = pi / 4

    return Args(current=current, length=length, angle=angle, induction=induction)


def test_basic_force(test_args: Args) -> None:
    result = ampere_force_from_current_and_induction.calculate_force(test_args.current, test_args.length, test_args.angle,
        test_args.induction)
    assert_equal(result, 3 * units.newton)


def test_bad_current(test_args: Args) -> None:
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ampere_force_from_current_and_induction.calculate_force(current, test_args.length, test_args.angle, test_args.induction)
    with raises(TypeError):
        ampere_force_from_current_and_induction.calculate_force(100, test_args.length, test_args.angle, test_args.induction)


def test_bad_length(test_args: Args) -> None:
    length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ampere_force_from_current_and_induction.calculate_force(test_args.current, length, test_args.angle, test_args.induction)
    with raises(TypeError):
        ampere_force_from_current_and_induction.calculate_force(test_args.current, 100, test_args.angle, test_args.induction)


def test_bad_angle(test_args: Args) -> None:
    angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ampere_force_from_current_and_induction.calculate_force(test_args.current, test_args.length, angle, test_args.induction)
    with raises(AttributeError):
        ampere_force_from_current_and_induction.calculate_force(test_args.current, test_args.length, True, test_args.induction)


def test_bad_induction(test_args: Args) -> None:
    induction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ampere_force_from_current_and_induction.calculate_force(test_args.current, test_args.length, test_args.angle, induction)
    with raises(TypeError):
        ampere_force_from_current_and_induction.calculate_force(test_args.current, test_args.length, test_args.angle, 100)
