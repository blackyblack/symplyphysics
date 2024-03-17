from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from sympy.physics.units import speed_of_light
from symplyphysics.laws.astronomy.relativistic import relative_speed_of_rocket_depends_on_mass_and_impulse as speed_law

# Description
## Let the specific impulse be equal to the speed of light. Then, with an initial mass of 72200 tons and
## a final mass of 200 tons, the speed will be 299787857 meters per second.
## https://scask.ru/d_book_msp.php?id=164

Args = namedtuple("Args", ["specific_impulse", "initial_mass", "final_mass"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # specific_impulse = Quantity(270000 * units.meter / units.second)
    specific_impulse = speed_of_light
    initial_mass = Quantity(72200 * units.tonne)
    final_mass = Quantity(200 * units.tonne)

    return Args(specific_impulse=specific_impulse, initial_mass=initial_mass, final_mass=final_mass)


def test_basic_speed(test_args: Args) -> None:
    result = speed_law.calculate_speed(test_args.specific_impulse, test_args.initial_mass,
        test_args.final_mass)
    assert_equal(result, 299787857 * units.meter / units.second)


def test_bad_specific_impulse(test_args: Args) -> None:
    specific_impulse = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(specific_impulse, test_args.initial_mass, test_args.final_mass)
    with raises(TypeError):
        speed_law.calculate_speed(100, test_args.initial_mass, test_args.final_mass)


def test_bad_mass(test_args: Args) -> None:
    bad_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(test_args.specific_impulse, bad_mass, test_args.final_mass)
    with raises(TypeError):
        speed_law.calculate_speed(test_args.specific_impulse, 100, test_args.final_mass)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(test_args.specific_impulse, test_args.initial_mass, bad_mass)
    with raises(TypeError):
        speed_law.calculate_speed(test_args.specific_impulse, test_args.initial_mass, 100)
