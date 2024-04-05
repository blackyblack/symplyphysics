from collections import namedtuple
from sympy import I
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.electricity import capacitor_impedance_from_capacitive_reactance as capacitor_impedance_law

# Description
## Assert we have a capacitor with capacitive reactance 1 Ohm.
## Then the capacitive impedance will be equal to -1 * I Ohm.

Args = namedtuple("Args", ["c"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c = Quantity(1 * units.ohm)
    return Args(c=c)


def test_basic_impedance(test_args: Args) -> None:
    result = capacitor_impedance_law.calculate_impedance(test_args.c)
    assert_equal(result, -1 * I * units.ohm)


def test_bad_reactance() -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitor_impedance_law.calculate_impedance(cb)
    with raises(TypeError):
        capacitor_impedance_law.calculate_impedance(100)
