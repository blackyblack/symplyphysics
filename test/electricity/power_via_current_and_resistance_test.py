from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import power_via_current_and_resistance as joule_lenz_law

# Description
## Assert we are having 1Amp flowing through 2-Ohm resistor.
## According to Joule-Lenz law we should have amount of heat dissipated on this resistor equals to 1^2 * 2 = 2 Watts.

Args = namedtuple("Args", ["C", "R"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = Quantity(1 * units.ampere)
    R = Quantity(2 * units.ohm)
    return Args(C=C, R=R)


def test_basic_power(test_args: Args) -> None:
    result = joule_lenz_law.calculate_heat_power(test_args.C, test_args.R)
    assert_equal(result, 2 * units.watt)


def test_bad_current(test_args: Args) -> None:
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_heat_power(Ib, test_args.R)
    with raises(TypeError):
        joule_lenz_law.calculate_heat_power(100, test_args.R)


def test_bad_resistance(test_args: Args) -> None:
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_heat_power(test_args.C, Rb)
    with raises(TypeError):
        joule_lenz_law.calculate_heat_power(test_args.C, 100)
