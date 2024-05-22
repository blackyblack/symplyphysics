from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits import internal_resistance_of_vacuum_diode as resistance_law

# Description
## The diode constant is 0.685 [milliampere / volt^(3 / 2)]. The voltage between the anode and the cathode is 17 volt.
## Then the internal resistance of the diode is 236 ohm.

Args = namedtuple("Args", ["diode_constant", "voltage"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    diode_constant = Quantity(0.685 * prefixes.milli * units.ampere / units.volt**(3 / 2))
    voltage = Quantity(17 * units.volt)

    return Args(diode_constant=diode_constant,
        voltage=voltage)


def test_basic_internal_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_internal_resistance(
        test_args.diode_constant, test_args.voltage)
    assert_equal(result, 236 * units.ohm)


def test_bad_diode_constant(test_args: Args) -> None:
    diode_constant = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_internal_resistance(diode_constant,
            test_args.voltage)
    with raises(TypeError):
        resistance_law.calculate_internal_resistance(100,
            test_args.voltage)


def test_bad_voltage(test_args: Args) -> None:
    voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_internal_resistance(
            test_args.diode_constant, voltage)
    with raises(TypeError):
        resistance_law.calculate_internal_resistance(
            test_args.diode_constant, 100)
