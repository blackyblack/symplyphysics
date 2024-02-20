from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.thermodynamics import gas_pressure_from_temperature as pressure_law

# Description
## With a pressure coefficient of 0.003661 [1 / kelvin], a temperature of 293.15 kelvin and a standard pressure of 101325 pascal,
## the pressure change will be 108744 pascal.
## https://www.indigomath.ru/formuly-po-fizike/molekuljarnaja-kinetika/temperaturnaja-zavisimost-davlenija-gaza.html
## https://www.fxyz.ru/%D1%84%D0%BE%D1%80%D0%BC%D1%83%D0%BB%D1%8B_%D0%BF%D0%BE_%D1%84%D0%B8%D0%B7%D0%B8%D0%BA%D0%B5/%D1%82%D0%B5%D1%80%D0%BC%D0%BE%D0%B4%D0%B8%D0%BD%D0%B0%D0%BC%D0%B8%D0%BA%D0%B0_%D1%82%D0%B5%D0%BE%D1%80%D0%B8%D1%8F_%D1%82%D0%B5%D0%BF%D0%BB%D0%BE%D1%82%D1%8B/%D1%82%D0%B5%D0%BF%D0%BB%D0%BE%D0%B2%D0%BE%D0%B5_%D1%80%D0%B0%D1%81%D1%88%D0%B8%D1%80%D0%B5%D0%BD%D0%B8%D0%B5/%D1%80%D0%B0%D1%81%D1%88%D0%B8%D1%80%D0%B5%D0%BD%D0%B8%D0%B5_%D0%B3%D0%B0%D0%B7%D0%B0/%D1%82%D0%B5%D1%80%D0%BC%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9_%D0%BA%D0%BE%D1%8D%D1%84%D1%84%D0%B8%D1%86%D0%B8%D0%B5%D0%BD%D1%82_%D0%B4%D0%B0%D0%B2%D0%BB%D0%B5%D0%BD%D0%B8%D1%8F/

Args = namedtuple("Args", ["standard_pressure", "thermal_coefficient", "temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    standard_pressure = Quantity(101325 * units.pascal)
    thermal_coefficient = Quantity(0.003661 * (1 / units.kelvin))
    temperature = Quantity(293.15 * units.kelvin)

    return Args(standard_pressure=standard_pressure,
        thermal_coefficient=thermal_coefficient,
        temperature=temperature)


def test_basic_pressure_change(test_args: Args) -> None:
    result = pressure_law.calculate_pressure_change(test_args.standard_pressure,
        test_args.thermal_coefficient, test_args.temperature)
    assert_equal(result, 108744 * units.pascal, tolerance=0.01)


def test_bad_pressure(test_args: Args) -> None:
    standard_pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_change(standard_pressure, test_args.thermal_coefficient,
            test_args.temperature)
    with raises(TypeError):
        pressure_law.calculate_pressure_change(100, test_args.thermal_coefficient,
            test_args.temperature)


def test_bad_thermal_coefficient(test_args: Args) -> None:
    thermal_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_change(test_args.standard_pressure, thermal_coefficient,
            test_args.temperature)
    with raises(TypeError):
        pressure_law.calculate_pressure_change(test_args.standard_pressure, 100,
            test_args.temperature)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_change(test_args.standard_pressure,
            test_args.thermal_coefficient, temperature)
    with raises(TypeError):
        pressure_law.calculate_pressure_change(test_args.standard_pressure,
            test_args.thermal_coefficient, 100)
