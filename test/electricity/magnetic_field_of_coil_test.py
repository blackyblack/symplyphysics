from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import magnetic_field_of_coil as induction_law

# Description
## With the number of turns equal to 11, a current of 1 ampere, a solenoid length of 0.5 meter
## the magnetic induction will be 2.76e-5 tesla.
## https://physics.icalculator.com/magnetic-field-inside-a-solenoid-calculator.html

Args = namedtuple("Args", ["current", "length", "number_turns"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    current = Quantity(1 * units.ampere)
    length = Quantity(0.5 * units.meter)
    number_turns = 11
    return Args(current=current, length=length, number_turns=number_turns)


def test_basic_induction(test_args: Args) -> None:
    result = induction_law.calculate_induction(test_args.current, test_args.length,
        test_args.number_turns)
    assert_equal(result, 2.76e-5 * units.tesla, relative_tolerance=2e-3)


def test_bad_current(test_args: Args) -> None:
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(current, test_args.length, test_args.number_turns)
    with raises(TypeError):
        induction_law.calculate_induction(100, test_args.length, test_args.number_turns)


def test_bad_length(test_args: Args) -> None:
    length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.current, length, test_args.number_turns)
    with raises(TypeError):
        induction_law.calculate_induction(test_args.current, 100, test_args.number_turns)


def test_bad_number_turns(test_args: Args) -> None:
    number_turns = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.current, test_args.length, number_turns)
    with raises(ValueError):
        induction_law.calculate_induction(test_args.current, test_args.length, -1)
