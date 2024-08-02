from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.laws.thermodynamics import internal_energy_of_ideal_gas_via_temperature as inner_energy_law

# Description
## With help of calculator on https://www.calculatoratoz.com/ru/molar-internal-energy-of-ideal-gas-calculator/Calc-1705 I calculated inner energy of 1 mole (4 gram) of Helium with it's molar mass 4 gramm/mole at 20 Celsius degrees.
## It should be 3656 Joules of energy.

Args = namedtuple("Args", ["m", "T", "M"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(4 * units.gram)
    normal_temperature = Celsius(20)
    T = to_kelvin_quantity(normal_temperature)
    M = Quantity(4 * units.gram / units.mole)
    return Args(m=m, T=T, M=M)


def test_basic_energy(test_args: Args) -> None:
    result = inner_energy_law.calculate_inner_energy(test_args.m, test_args.T, test_args.M)
    assert_equal(result, 3656 * units.joule)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inner_energy_law.calculate_inner_energy(mb, test_args.T, test_args.M)
    with raises(TypeError):
        inner_energy_law.calculate_inner_energy(100, test_args.T, test_args.M)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inner_energy_law.calculate_inner_energy(test_args.m, tb, test_args.M)
    with raises(TypeError):
        inner_energy_law.calculate_inner_energy(test_args.m, 100, test_args.M)


def test_bad_mole_mass(test_args: Args) -> None:
    Mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inner_energy_law.calculate_inner_energy(test_args.m, test_args.T, Mb)
    with raises(TypeError):
        inner_energy_law.calculate_inner_energy(test_args.m, test_args.T, 100)
