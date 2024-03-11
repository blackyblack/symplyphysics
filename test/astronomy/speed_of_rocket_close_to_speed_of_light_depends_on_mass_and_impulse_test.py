from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import speed_of_rocket_close_to_speed_of_light_depends_on_mass_and_impulse as speed_law

# Description
## With the initial mass of the rocket equal to 100 tons, the final mass of the rocket equal to 20 tons,
## and a specific impulse equal to 3000 meters per second, the rocket speed will be 4828.31 meters per second.

Args = namedtuple("Args", ["specific_impulse", "initial_mass", "final_mass"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    specific_impulse = Quantity(3000 * units.meter / units.second)
    initial_mass = Quantity(100 * units.tonne)
    final_mass = Quantity(20 * units.tonne)

    return Args(specific_impulse=specific_impulse,
        initial_mass=initial_mass,
        final_mass=final_mass)


def test_basic_speed(test_args: Args) -> None:
    result = speed_law.calculate_speed(test_args.specific_impulse,
        test_args.initial_mass, test_args.final_mass)
    assert_equal(result, 4828.31 * units.meter / units.second)


def test_bad_specific_impulse(test_args: Args) -> None:
    specific_impulse = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(specific_impulse, test_args.initial_mass,
            test_args.final_mass)
    with raises(TypeError):
        speed_law.calculate_speed(100, test_args.initial_mass, test_args.final_mass)


def test_bad_mass(test_args: Args) -> None:
    bad_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(test_args.specific_impulse, bad_mass,
            test_args.final_mass)
    with raises(TypeError):
        speed_law.calculate_speed(test_args.specific_impulse, 100,
            test_args.final_mass)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(test_args.specific_impulse,
            test_args.initial_mass, bad_mass)
    with raises(TypeError):
        speed_law.calculate_speed(test_args.specific_impulse,
            test_args.initial_mass, 100)
