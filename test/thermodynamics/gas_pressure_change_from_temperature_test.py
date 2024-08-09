from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, Quantity, errors, prefixes
from symplyphysics.laws.thermodynamics import gas_pressure_change_from_temperature as pressure_law

# Description
## With a temperature of 293.15 kelvin and a standard pressure of 101325 pascal,
## the pressure change will be 7.42 kPa.
## https://www.indigomath.ru/formuly-po-fizike/molekuljarnaja-kinetika/temperaturnaja-zavisimost-davlenija-gaza.html

Args = namedtuple("Args", ["standard_pressure", "temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    standard_pressure = Quantity(101325 * units.pascal)
    temperature = Quantity(293.15 * units.kelvin)
    return Args(standard_pressure=standard_pressure, temperature=temperature)


def test_basic_pressure_change(test_args: Args) -> None:
    result = pressure_law.calculate_pressure_change(test_args.standard_pressure,
        test_args.temperature)
    assert_equal(result, 7.42 * prefixes.kilo * units.pascal)


def test_bad_pressure(test_args: Args) -> None:
    standard_pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_change(standard_pressure, test_args.temperature)
    with raises(TypeError):
        pressure_law.calculate_pressure_change(100, test_args.temperature)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_change(test_args.standard_pressure, temperature)
    with raises(TypeError):
        pressure_law.calculate_pressure_change(test_args.standard_pressure, 100)
