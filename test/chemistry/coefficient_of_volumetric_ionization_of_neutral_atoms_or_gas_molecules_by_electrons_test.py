from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.chemistry import coefficient_of_volumetric_ionization_of_neutral_atoms_or_gas_molecules_by_electrons as coefficient_law

## The first constant of gas is equal to 0.099 [1 / centimeter / pascal]. The second constant of gas is 1.84 [volt / centimeter / pascal].
## The pressure is 293 pascal. The electric intensity is 659 [volt / centimeter].
## Then the coefficient of volumetric ionization will be 12.8 [1 / centimeter].

Args = namedtuple("Args", [
    "first_constant_of_gas",
    "second_constant_of_gas",
    "pressure",
    "electric_intensity",
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    first_constant_of_gas = Quantity(0.099 / units.centimeter / units.pascal)
    second_constant_of_gas = Quantity(1.84 * units.volt / units.centimeter / units.pascal)
    pressure = Quantity(293 * units.pascal)
    electric_intensity = Quantity(659 * units.volt / units.centimeter)

    return Args(first_constant_of_gas=first_constant_of_gas,
        second_constant_of_gas=second_constant_of_gas,
        pressure=pressure,
        electric_intensity=electric_intensity)


def test_basic_coefficient_of_volumetric_ionization(test_args: Args) -> None:
    result = coefficient_law.calculate_coefficient_of_volumetric_ionization(
        test_args.first_constant_of_gas, test_args.second_constant_of_gas, test_args.pressure,
        test_args.electric_intensity)
    assert_equal(result, 12.8 * (1 / units.centimeter))


def test_bad_first_constant_of_gas(test_args: Args) -> None:
    first_constant_of_gas = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_of_volumetric_ionization(
            first_constant_of_gas, test_args.second_constant_of_gas, test_args.pressure,
            test_args.electric_intensity)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_of_volumetric_ionization(
            100, test_args.second_constant_of_gas, test_args.pressure, test_args.electric_intensity)


def test_bad_second_constant_of_gas(test_args: Args) -> None:
    second_constant_of_gas = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_of_volumetric_ionization(
            test_args.first_constant_of_gas, second_constant_of_gas, test_args.pressure,
            test_args.electric_intensity)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_of_volumetric_ionization(
            test_args.first_constant_of_gas, 100, test_args.pressure, test_args.electric_intensity)


def test_bad_pressure(test_args: Args) -> None:
    pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_of_volumetric_ionization(
            test_args.first_constant_of_gas, test_args.second_constant_of_gas, pressure,
            test_args.electric_intensity)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_of_volumetric_ionization(
            test_args.first_constant_of_gas, test_args.second_constant_of_gas, 100,
            test_args.electric_intensity)


def test_bad_electric_intensity(test_args: Args) -> None:
    electric_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_of_volumetric_ionization(
            test_args.first_constant_of_gas, test_args.second_constant_of_gas, test_args.pressure,
            electric_intensity)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_of_volumetric_ionization(
            test_args.first_constant_of_gas, test_args.second_constant_of_gas, test_args.pressure,
            100)
