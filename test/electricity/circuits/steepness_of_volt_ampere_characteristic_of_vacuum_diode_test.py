from collections import namedtuple
from pytest import fixture, raises
from sympy import Rational
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits import steepness_of_volt_ampere_characteristic_of_vacuum_diode as steepness_law

# Description
## The diode constant is 0.685 [milliampere / volt^(3 / 2)]. The voltage between the anode and the cathode is 17 volt.
## Then the steepness of the volt-ampere characteristic of the diode is 4.24 milliampere per volt.

Args = namedtuple("Args", ["diode_constant", "voltage"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    diode_constant = Quantity(0.685 * prefixes.milli * units.ampere / units.volt**Rational(3, 2))
    voltage = Quantity(17 * units.volt)

    return Args(diode_constant=diode_constant,
        voltage=voltage)


def test_basic_steepness(test_args: Args) -> None:
    result = steepness_law.calculate_steepness(
        test_args.diode_constant, test_args.voltage)
    assert_equal(result, 4.24 * prefixes.milli * units.ampere / units.volt)


def test_bad_diode_constant(test_args: Args) -> None:
    diode_constant = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        steepness_law.calculate_steepness(diode_constant,
            test_args.voltage)
    with raises(TypeError):
        steepness_law.calculate_steepness(100,
            test_args.voltage)


def test_bad_voltage(test_args: Args) -> None:
    voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        steepness_law.calculate_steepness(
            test_args.diode_constant, voltage)
    with raises(TypeError):
        steepness_law.calculate_steepness(
            test_args.diode_constant, 100)
