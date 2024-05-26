from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.chemistry import breakdown_voltage_of_gas_discharge as voltage_law

## The first constant of gas is equal to 0.099 [1 / centimeter / pascal]. The second constant of gas is 1.84 [volt / centimeter / pascal].
## The pressure is 293 pascal. The distance between electrodes is equal to 1 centimeter. 
## The secondary emission factor of the ion-electronic type is 0.022.
## Then the breakdown voltage of the gas discharge will be 266.57 volt.

Args = namedtuple("Args", [
    "first_constant_of_gas", "second_constant_of_gas", "pressure", "distance_between_electrodes", "secondary_emission_factor"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    first_constant_of_gas = Quantity(0.099 / units.centimeter / units.pascal)
    second_constant_of_gas = Quantity(1.84 * units.volt / units.centimeter / units.pascal)
    pressure = Quantity(293 * units.pascal)
    distance_between_electrodes = Quantity(1 * units.centimeter)
    secondary_emission_factor = 0.022

    return Args(first_constant_of_gas=first_constant_of_gas,
        second_constant_of_gas=second_constant_of_gas,
        pressure=pressure,
        distance_between_electrodes=distance_between_electrodes,
        secondary_emission_factor=secondary_emission_factor)


def test_basic_breakdown_voltage(test_args: Args) -> None:
    result = voltage_law.calculate_breakdown_voltage(
        test_args.first_constant_of_gas, test_args.second_constant_of_gas, test_args.pressure,
        test_args.distance_between_electrodes, test_args.secondary_emission_factor,)
    assert_equal(result, 266.57 * units.volt)


def test_bad_first_constant_of_gas(test_args: Args) -> None:
    first_constant_of_gas = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_breakdown_voltage(first_constant_of_gas,
            test_args.second_constant_of_gas, test_args.pressure,
            test_args.distance_between_electrodes, test_args.secondary_emission_factor)
    with raises(TypeError):
        voltage_law.calculate_breakdown_voltage(100,
            test_args.second_constant_of_gas, test_args.pressure,
            test_args.distance_between_electrodes, test_args.secondary_emission_factor)


def test_bad_second_constant_of_gas(test_args: Args) -> None:
    second_constant_of_gas = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_breakdown_voltage(test_args.first_constant_of_gas,
            second_constant_of_gas, test_args.pressure, test_args.distance_between_electrodes, test_args.secondary_emission_factor)
    with raises(TypeError):
        voltage_law.calculate_breakdown_voltage(test_args.first_constant_of_gas,
            100, test_args.pressure, test_args.distance_between_electrodes, test_args.secondary_emission_factor)


def test_bad_pressure(test_args: Args) -> None:
    pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_breakdown_voltage(test_args.first_constant_of_gas,
            test_args.second_constant_of_gas, pressure, test_args.distance_between_electrodes, test_args.secondary_emission_factor)
    with raises(TypeError):
        voltage_law.calculate_breakdown_voltage(test_args.first_constant_of_gas,
            test_args.second_constant_of_gas, 100, test_args.distance_between_electrodes, test_args.secondary_emission_factor)


def test_bad_distance_between_electrodes(test_args: Args) -> None:
    distance_between_electrodes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_breakdown_voltage(test_args.first_constant_of_gas,
            test_args.second_constant_of_gas, test_args.pressure, distance_between_electrodes, test_args.secondary_emission_factor)
    with raises(TypeError):
        voltage_law.calculate_breakdown_voltage(test_args.first_constant_of_gas,
            test_args.second_constant_of_gas, test_args.pressure, 100, test_args.secondary_emission_factor)


def test_bad_secondary_emission_factor(test_args: Args) -> None:
    secondary_emission_factor = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_breakdown_voltage(test_args.first_constant_of_gas,
            test_args.second_constant_of_gas, test_args.pressure, test_args.distance_between_electrodes, secondary_emission_factor)
