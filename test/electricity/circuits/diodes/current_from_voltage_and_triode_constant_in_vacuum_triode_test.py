from collections import namedtuple
from pytest import fixture, raises
from sympy import Rational
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits.diodes import current_from_voltage_and_triode_constant_in_vacuum_triode as current_law

# Description
## The constant of the vacuum triode is 0.685 [milliampere / volt^(3 / 2)]. The voltage between the anode and the cathode
## is 17 volt. The voltage triode gain is 3. The grid voltage is -2 volt. Then the anode current will be equal to 25 milliampere.

Args = namedtuple("Args",
    ["triode_constant", "anode_voltage", "voltage_triode_gain", "grid_voltage"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    triode_constant = Quantity(0.685 * prefixes.milli * units.ampere / units.volt**Rational(3, 2))
    anode_voltage = Quantity(17 * units.volt)
    voltage_triode_gain = 3
    grid_voltage = Quantity(-2 * units.volt)

    return Args(triode_constant=triode_constant,
        anode_voltage=anode_voltage,
        voltage_triode_gain=voltage_triode_gain,
        grid_voltage=grid_voltage)


def test_basic_current(test_args: Args) -> None:
    result = current_law.calculate_current(test_args.triode_constant, test_args.anode_voltage,
        test_args.voltage_triode_gain, test_args.grid_voltage)
    assert_equal(result, 25 * prefixes.milli * units.ampere)


def test_bad_triode_constant(test_args: Args) -> None:
    triode_constant = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(triode_constant, test_args.anode_voltage,
            test_args.voltage_triode_gain, test_args.grid_voltage)
    with raises(TypeError):
        current_law.calculate_current(100, test_args.anode_voltage, test_args.voltage_triode_gain,
            test_args.grid_voltage)


def test_bad_anode_voltage(test_args: Args) -> None:
    anode_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.triode_constant, anode_voltage,
            test_args.voltage_triode_gain, test_args.grid_voltage)
    with raises(TypeError):
        current_law.calculate_current(test_args.triode_constant, 100, test_args.voltage_triode_gain,
            test_args.grid_voltage)


def test_bad_voltage_triode_gain(test_args: Args) -> None:
    voltage_triode_gain = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.triode_constant, test_args.anode_voltage,
            voltage_triode_gain, test_args.grid_voltage)


def test_bad_grid_voltage(test_args: Args) -> None:
    grid_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.triode_constant, test_args.anode_voltage,
            test_args.voltage_triode_gain, grid_voltage)
    with raises(TypeError):
        current_law.calculate_current(test_args.triode_constant, test_args.anode_voltage,
            test_args.voltage_triode_gain, 100)
