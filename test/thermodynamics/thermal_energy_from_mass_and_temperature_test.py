from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as amount_energy

# How much energy does it take to heat some volume of water weighing 0.5 kilograms то 50 Celsius degree?
# Specific heat capacity of water is 4200 J/kg*K , ignore losses.

Args = namedtuple("Args", ["C", "m", "t1", "t2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = Quantity(4.2 * prefixes.kilo * units.joule / (units.kilogram * units.kelvin))
    m = Quantity(0.5 * units.kilogram)
    initial_temperature = Celsius(0)
    t1 = to_kelvin_quantity(initial_temperature)
    final_temperature = Celsius(50)
    t2 = to_kelvin_quantity(final_temperature)
    return Args(C=C, m=m, t1=t1, t2=t2)


def test_basic_amount(test_args: Args) -> None:
    result = amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2,
        test_args.t1)
    assert_equal(result, 105000.1 * units.joule)


def test_bad_specific_heat(test_args: Args) -> None:
    Cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(Cb, test_args.m, test_args.t2, test_args.t1)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(100, test_args.m, test_args.t2, test_args.t1)


def test_bad_body_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, mb, test_args.t2, test_args.t1)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, 100, test_args.t2, test_args.t1)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, tb, test_args.t1)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, 100, test_args.t1)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2, tb)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2, 100)
