from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import magnetic_induction_of_solenoid as induction_law

# Description
## With the number of turns equal to 11, a current of 1 ampere, a solenoid length of 0.5 meter
## and a magnetic permeability of 120, the magnetic induction will be 3.3174e-3 tesla.
## https://physics.icalculator.com/magnetic-field-inside-a-solenoid-calculator.html


@fixture(name="test_args")
def test_args_fixture():
    current = Quantity(1 * units.ampere)
    length = Quantity(0.5 * units.meter)
    number_turns = 11
    relative_permeability = 120

    Args = namedtuple("Args", ["current", "length", "number_turns", "relative_permeability"])
    return Args(current=current,
        length=length,
        relative_permeability=relative_permeability,
        number_turns=number_turns)


def test_basic_induction(test_args):
    result = induction_law.calculate_induction(test_args.current, test_args.length,
        test_args.relative_permeability, test_args.number_turns)
    assert_equal(result, 3.3174e-3 * units.tesla)


def test_bad_current(test_args):
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(current, test_args.length,
            test_args.relative_permeability, test_args.number_turns)
    with raises(TypeError):
        induction_law.calculate_induction(100, test_args.length, test_args.relative_permeability,
            test_args.number_turns)


def test_bad_length(test_args):
    length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.current, length,
            test_args.relative_permeability, test_args.number_turns)
    with raises(TypeError):
        induction_law.calculate_induction(test_args.current, 100, test_args.relative_permeability,
            test_args.number_turns)


def test_bad_relative_permeability(test_args):
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.current, test_args.length,
            relative_permeability, test_args.number_turns)


def test_bad_number_turns(test_args):
    number_turns = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.current, test_args.length,
            test_args.relative_permeability, number_turns)
    with raises(ValueError):
        induction_law.calculate_induction(test_args.current, test_args.length,
            test_args.relative_permeability, -1)
