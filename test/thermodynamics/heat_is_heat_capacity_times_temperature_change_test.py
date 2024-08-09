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
from symplyphysics.laws.thermodynamics import heat_is_heat_capacity_times_temperature_change as thermal_energy

# How much energy does it take to heat some volume of water tо 50 degree Celsius?
# Heat capacity of that amount of water is 2.1 kJ/K, ignore losses.

Args = namedtuple("Args", ["C", "dt"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = Quantity(2.1 * prefixes.kilo * units.joule / units.kelvin)
    initial_temperature = Celsius(0)
    t1 = to_kelvin_quantity(initial_temperature)
    final_temperature = Celsius(50)
    t2 = to_kelvin_quantity(final_temperature)
    dt = Quantity(t2 - t1)
    return Args(C=C, dt=dt)


def test_basic_amount(test_args: Args) -> None:
    result = thermal_energy.calculate_amount_energy(test_args.C, test_args.dt)
    assert_equal(result, 105000.1 * units.joule)


def test_bad_heat_capacity(test_args: Args) -> None:
    Cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thermal_energy.calculate_amount_energy(Cb, test_args.dt)
    with raises(TypeError):
        thermal_energy.calculate_amount_energy(100, test_args.dt)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thermal_energy.calculate_amount_energy(test_args.C, tb)
    with raises(TypeError):
        thermal_energy.calculate_amount_energy(test_args.C, 100)
