from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import electron_distribution_function_in_gas_plasma_by_druyvestein as function_law

# Description
## The voltage between electrodes is 250 volt. The electron energy is equal to 110 electronvolt.
## Then the value of the electron distribution function is 8.32e-4.

Args = namedtuple("Args", ["voltage_between_electrodes", "electron_energy"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    voltage_between_electrodes = Quantity(250 * units.volt)
    electron_energy = Quantity(110 * units.electronvolt)

    return Args(voltage_between_electrodes=voltage_between_electrodes,
        electron_energy=electron_energy)


def test_basic_value_of_distribution_function(test_args: Args) -> None:
    result = function_law.calculate_value_of_distribution_function(
        test_args.voltage_between_electrodes, test_args.electron_energy)
    assert_equal(result, 8.32e-4)

def test_bad_voltage_between_electrodes(test_args: Args) -> None:
    voltage_between_electrodes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        function_law.calculate_value_of_distribution_function(voltage_between_electrodes,
            test_args.electron_energy)
    with raises(TypeError):
        function_law.calculate_value_of_distribution_function(100,
            test_args.electron_energy)


def test_bad_electron_energy(test_args: Args) -> None:
    electron_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        function_law.calculate_value_of_distribution_function(
            test_args.voltage_between_electrodes, electron_energy)
    with raises(TypeError):
        function_law.calculate_value_of_distribution_function(
            test_args.voltage_between_electrodes, 100)
