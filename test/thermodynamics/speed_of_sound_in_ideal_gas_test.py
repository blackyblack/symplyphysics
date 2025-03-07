from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.laws.thermodynamics import speed_of_sound_in_ideal_gas

# Description
# Input: Air, temperature=20°C, gamma=1.4, M=29 g/mol
# Comparing with the tabular value from Wikipedia
# It should be 343.21 m/s

Args = namedtuple("Args", ["t", "gamma", "M"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    gamma = 1.4
    t = to_kelvin_quantity(Celsius(20))
    M = Quantity(29 * units.gram / units.mole)
    return Args(t=t, gamma=gamma, M=M)


def test_speed_of_sound(test_args: Args) -> None:
    result = speed_of_sound_in_ideal_gas.calculate_speed_of_sound(test_args.t, test_args.gamma,
        test_args.M)
    assert_equal(result, 343.21 * units.meter / units.second)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_of_sound_in_ideal_gas.calculate_speed_of_sound(test_args.t, test_args.gamma, tb)
    with raises(TypeError):
        speed_of_sound_in_ideal_gas.calculate_speed_of_sound(test_args.t, test_args.gamma, 100)


def test_bad_mole_mass(test_args: Args) -> None:
    Mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_of_sound_in_ideal_gas.calculate_speed_of_sound(test_args.t, test_args.gamma, Mb)
    with raises(TypeError):
        speed_of_sound_in_ideal_gas.calculate_speed_of_sound(test_args.t, test_args.gamma, 100)


def test_bad_heat_capacity_ratio(test_args: Args) -> None:
    gamma = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_of_sound_in_ideal_gas.calculate_speed_of_sound(test_args.t, gamma, test_args.M)
