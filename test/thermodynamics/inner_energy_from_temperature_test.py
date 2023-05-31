from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.laws.thermodynamics import inner_energy_from_temperature as inner_energy_law

# Description
## With help of calculator on https://www.calculatoratoz.com/ru/molar-internal-energy-of-ideal-gas-calculator/Calc-1705 I calculated inner energy of 1 mole (4 gram) of Helium with it's molar mass 4 gramm/mole at 20 Celsius degrees.
## It should be 3656 Joules of energy.


@fixture
def test_args():
    m = Quantity(4 * units.gram)
    normal_temperature = Celsius(20)
    T = to_kelvin_quantity(normal_temperature)
    M = Quantity(4 * units.gram / units.mole)
    Args = namedtuple("Args", ["m", "T", "M"])
    return Args(m=m, T=T, M=M)


def test_basic_energy(test_args):
    result = inner_energy_law.calculate_inner_energy(test_args.m, test_args.T, test_args.M)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_pressure = convert_to(result, units.joule).subs(units.joule, 1).evalf(2)
    assert result_pressure == approx(3656, 0.01)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inner_energy_law.calculate_inner_energy(mb, test_args.T, test_args.M)
    with raises(TypeError):
        inner_energy_law.calculate_inner_energy(100, test_args.T, test_args.M)


def test_bad_temperature(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inner_energy_law.calculate_inner_energy(test_args.m, tb, test_args.M)
    with raises(TypeError):
        inner_energy_law.calculate_inner_energy(test_args.m, 100, test_args.M)


def test_bad_mole_mass(test_args):
    Mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inner_energy_law.calculate_inner_energy(test_args.m, test_args.T, Mb)
    with raises(TypeError):
        inner_energy_law.calculate_inner_energy(test_args.m, test_args.T, 100)
