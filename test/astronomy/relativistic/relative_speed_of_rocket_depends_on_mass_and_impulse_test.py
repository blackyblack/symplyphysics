from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy.relativistic import relative_speed_of_rocket_depends_on_mass_and_impulse as speed_law

# Description
## Let the exhaust velocity be equal to the speed of light. Then, with an initial mass of 72200 tons and
## a final mass of 200 tons, the speed will be 299787857 meters per second.
## https://scask.ru/d_book_msp.php?id=164

Args = namedtuple("Args", ["exhaust_velocity", "initial_mass", "final_mass"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    exhaust_velocity = speed_of_light
    initial_mass = Quantity(72200 * units.tonne)
    final_mass = Quantity(200 * units.tonne)

    return Args(exhaust_velocity=exhaust_velocity, initial_mass=initial_mass, final_mass=final_mass)


def test_basic_speed(test_args: Args) -> None:
    result = speed_law.calculate_speed(test_args.exhaust_velocity, test_args.initial_mass,
        test_args.final_mass)
    assert_equal(result, 299787857 * units.meter / units.second)


def test_bad_exhaust_velocity(test_args: Args) -> None:
    exhaust_velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(exhaust_velocity, test_args.initial_mass, test_args.final_mass)
    with raises(TypeError):
        speed_law.calculate_speed(100, test_args.initial_mass, test_args.final_mass)


def test_bad_mass(test_args: Args) -> None:
    bad_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(test_args.exhaust_velocity, bad_mass, test_args.final_mass)
    with raises(TypeError):
        speed_law.calculate_speed(test_args.exhaust_velocity, 100, test_args.final_mass)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(test_args.exhaust_velocity, test_args.initial_mass, bad_mass)
    with raises(TypeError):
        speed_law.calculate_speed(test_args.exhaust_velocity, test_args.initial_mass, 100)
