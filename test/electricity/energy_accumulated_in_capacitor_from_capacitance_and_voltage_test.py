from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity, prefixes)
from symplyphysics.laws.electricity import energy_accumulated_in_capacitor_from_capacitance_and_voltage as capacitor_law

# Description
## Assert we have 220 uF capacitor charged to 10 volts.
## According to law we should have amount of energy accumulated in this capacitor equals to 220 * 0.000001 * 10**2 / 2 = 0.011 Joules.

Args = namedtuple("Args", ["C", "V"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = Quantity(220 * prefixes.micro * units.farad)
    V = Quantity(10 * units.volt)
    return Args(C=C, V=V)


def test_basic_energy(test_args: Args) -> None:
    result = capacitor_law.calculate_accumulated_energy(test_args.C, test_args.V)
    assert_equal(result, 0.011 * units.joule)


def test_bad_capacitance(test_args: Args) -> None:
    Cb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitor_law.calculate_accumulated_energy(Cb, test_args.V)
    with raises(TypeError):
        capacitor_law.calculate_accumulated_energy(100, test_args.V)


def test_bad_voltage(test_args: Args) -> None:
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitor_law.calculate_accumulated_energy(test_args.C, Vb)
    with raises(TypeError):
        capacitor_law.calculate_accumulated_energy(test_args.C, 100)
