# Description
## Assume we have a capacitor with 6 Coulombs of charge in it and having 3 volts on it's terminals.
## According to the definition of Capacitance this capacitor should have C = 2 Farads.

from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_def

Args = namedtuple("Args", ["C", "Q", "U"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = Quantity(2 * units.farad)
    Q = Quantity(6 * units.coulomb)
    U = Quantity(3 * units.volts)
    return Args(C=C, Q=Q, U=U)


def test_basic_capacitance(test_args: Args) -> None:
    result = capacitance_def.calculate_capacitance(test_args.Q, test_args.U)
    assert_equal(result, 2 * units.farad)


def test_bad_charge(test_args: Args) -> None:
    Qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitance_def.calculate_capacitance(Qb, test_args.U)
    with raises(TypeError):
        capacitance_def.calculate_capacitance(100, test_args.U)


def test_bad_voltage(test_args: Args) -> None:
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitance_def.calculate_capacitance(test_args.Q, Vb)
    with raises(TypeError):
        capacitance_def.calculate_capacitance(test_args.Q, 100)
