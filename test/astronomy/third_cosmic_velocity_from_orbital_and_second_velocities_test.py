from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors,)
from symplyphysics.laws.astronomy import third_cosmic_velocity_from_orbital_and_second_velocities as velocity_law

# Description
## The earth's orbital velocity is 29.783 kilometers per second. The second cosmic velocity of the Earth is 11.183 kilometers
## per second. Then the third cosmic velocity of the Earth is 16.65 kilometers per second.
## https://ru.wikipedia.org/wiki/Третья_космическая_скорость

Args = namedtuple("Args", ["orbital_velocity", "second_velocity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    orbital_velocity = Quantity(29.783 * (units.kilometer / units.second))
    second_velocity = Quantity(11.183 * (units.kilometer / units.second))

    return Args(orbital_velocity=orbital_velocity, second_velocity=second_velocity)


def test_basic_third_velocity(test_args: Args) -> None:
    result = velocity_law.calculate_third_velocity(test_args.orbital_velocity, test_args.second_velocity)
    assert_equal(result, 16.65 * (units.kilometer / units.second))


def test_bad_velocities(test_args: Args) -> None:
    bad_velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_third_velocity(bad_velocity, test_args.second_velocity)
    with raises(TypeError):
        velocity_law.calculate_third_velocity(100, test_args.second_velocity)
    with raises(errors.UnitsError):
        velocity_law.calculate_third_velocity(test_args.orbital_velocity, bad_velocity)
    with raises(TypeError):
        velocity_law.calculate_third_velocity(test_args.orbital_velocity, 100)
