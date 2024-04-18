from collections import namedtuple
from sympy import I
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.electricity import coil_impedance_from_inductive_reactance as coil_impedance_law

# Description
## Assert we have a coil with inductive reactance 1 Ohm.
## Then the coil impedance will be equal to 1 * I Ohm.

Args = namedtuple("Args", ["l"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    l = Quantity(1 * units.ohm)
    return Args(l=l)


def test_basic_impedance(test_args: Args) -> None:
    result = coil_impedance_law.calculate_impedance(test_args.l)
    assert_equal(result, 1 * I * units.ohm)


def test_bad_reactance() -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coil_impedance_law.calculate_impedance(cb)
    with raises(TypeError):
        coil_impedance_law.calculate_impedance(100)
