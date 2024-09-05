from collections import namedtuple
from sympy import I
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.electricity.circuits import coil_impedance_via_inductance_and_frequency as coil_impedance_law

# Description
## Assert we have a coil with 5H inductivity. In 50Hz (314.159 rad/s) circuit impedance of this coil should be 1570.8 Ohm.
## (https://www.translatorscafe.com/unit-converter/ru-RU/calculator/inductor-impedance/)

Args = namedtuple("Args", ["inductivity", "frequency"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    inductivity = Quantity(5 * units.henry)
    frequency = Quantity(314.159 * units.radian / units.second)
    return Args(inductivity=inductivity, frequency=frequency)


def test_basic_impedance(test_args: Args) -> None:
    result = coil_impedance_law.calculate_impedance(test_args.inductivity, test_args.frequency)
    assert_equal(result, 1570.8 * I * units.ohm)


def test_bad_inductivity(test_args: Args) -> None:
    ib = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coil_impedance_law.calculate_impedance(ib, test_args.frequency)
    with raises(TypeError):
        coil_impedance_law.calculate_impedance(100, test_args.frequency)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coil_impedance_law.calculate_impedance(test_args.inductivity, wb)
    with raises(TypeError):
        coil_impedance_law.calculate_impedance(test_args.inductivity, 100)
