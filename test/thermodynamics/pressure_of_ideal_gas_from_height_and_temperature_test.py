from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.thermodynamics import pressure_of_ideal_gas_from_height_and_temperature as pressure_of_ideal_gas

## Source of numbers: https://tsput.ru/res/fizika/1/ZAD_MKT/z_mkt_18_02.htm

Args = namedtuple("Args", ["initial_pressure", "molecular_mass", "height_change", "temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    initial_pressure = Quantity(1e5 * units.pascal)
    atomic_weight = Quantity(0.029 * units.kilogram / units.mole)
    molecular_mass = Quantity(atomic_weight / units.avogadro_constant)
    final_height = Quantity(6000 * units.meter)
    initial_height = Quantity(1000 * units.meter)
    height_change = Quantity(final_height - initial_height)
    temperature = Quantity(293 * units.kelvin)
    return Args(initial_pressure=initial_pressure,
        molecular_mass=molecular_mass,
        height_change=height_change,
        temperature=temperature)


def test_basic_law(test_args: Args) -> None:
    result = pressure_of_ideal_gas.calculate_final_pressure(test_args.initial_pressure,
        test_args.molecular_mass, test_args.height_change, test_args.temperature)
    assert_equal(result, 5.6e4 * units.pascal, relative_tolerance=0.01)


def test_bad_pressure(test_args: Args) -> None:
    bad_pressure = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        pressure_of_ideal_gas.calculate_final_pressure(bad_pressure, test_args.molecular_mass,
            test_args.height_change, test_args.temperature)
    with raises(TypeError):
        pressure_of_ideal_gas.calculate_final_pressure(100, test_args.molecular_mass,
            test_args.height_change, test_args.temperature)


def test_bad_molecular_mass(test_args: Args) -> None:
    bad_molecular_mass = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        pressure_of_ideal_gas.calculate_final_pressure(test_args.initial_pressure,
            bad_molecular_mass, test_args.height_change, test_args.temperature)
    with raises(TypeError):
        pressure_of_ideal_gas.calculate_final_pressure(test_args.initial_pressure, 100,
            test_args.height_change, test_args.temperature)


def test_bad_height(test_args: Args) -> None:
    bad_height = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_of_ideal_gas.calculate_final_pressure(test_args.initial_pressure,
            test_args.molecular_mass, bad_height, test_args.temperature)
    with raises(TypeError):
        pressure_of_ideal_gas.calculate_final_pressure(test_args.initial_pressure,
            test_args.molecular_mass, 100, test_args.temperature)


def test_bad_temperature(test_args: Args) -> None:
    bad_temperature = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        pressure_of_ideal_gas.calculate_final_pressure(test_args.initial_pressure,
            test_args.molecular_mass, test_args.height_change, bad_temperature)
    with raises(TypeError):
        pressure_of_ideal_gas.calculate_final_pressure(test_args.initial_pressure,
            test_args.molecular_mass, test_args.height_change, 100)
