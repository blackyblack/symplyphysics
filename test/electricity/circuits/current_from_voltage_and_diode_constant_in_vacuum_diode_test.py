from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits import current_from_voltage_and_diode_constant_in_vacuum_diode as current_law

# Description
## The constant of the vacuum diode is 0.685 [milliampere / volt^(3 / 2)]. The voltage between the anode and the cathode
## is 17 volt. Then the anode current will be equal to 48 milliampere.

Args = namedtuple("Args", ["diode_constant", "voltage"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    diode_constant = Quantity(0.685 * prefixes.milli * units.ampere / units.volt**(3 / 2))
    voltage = Quantity(17 * units.volt)

    return Args(diode_constant=diode_constant,
        voltage=voltage)


def test_basic_current(test_args: Args) -> None:
    result = current_law.calculate_current(
        test_args.diode_constant, test_args.voltage)
    assert_equal(result, 48 * prefixes.milli * units.ampere)


def test_bad_diode_constant(test_args: Args) -> None:
    diode_constant = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(diode_constant,
            test_args.voltage)
    with raises(TypeError):
        current_law.calculate_current(100,
            test_args.voltage)


def test_bad_voltage(test_args: Args) -> None:
    voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(
            test_args.diode_constant, voltage)
    with raises(TypeError):
        current_law.calculate_current(
            test_args.diode_constant, 100)
