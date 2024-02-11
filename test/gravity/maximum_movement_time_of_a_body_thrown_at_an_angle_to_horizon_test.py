from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.gravity import maximum_movement_time_of_a_body_thrown_at_an_angle_to_horizon as time_law

# Description
## Let the initial velocity be 3 meter per second, the angle of the throw is 45 degree (pi / 4 radian),
## and the acceleration of gravity is 9.8 [meter / second^2]. Then the travel time will be 0.433 second.
## https://www.indigomath.ru//raschety/yJTxXc.html

Args = namedtuple("Args", ["initial_velocity", "angle"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    initial_velocity = Quantity(3 * (units.meter / units.second))
    angle = pi / 4
    return Args(initial_velocity=initial_velocity, angle=angle)


def test_basic_movement_time(test_args: Args) -> None:
    result = time_law.calculate_movement_time(test_args.initial_velocity, test_args.angle)
    assert_equal(result, 0.433 * units.second)


def test_bad_initial_velocity(test_args: Args) -> None:
    initial_velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_movement_time(initial_velocity, test_args.angle)
    with raises(TypeError):
        time_law.calculate_movement_time(100, test_args.angle)


def test_bad_angle(test_args: Args) -> None:
    angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_movement_time(test_args.initial_velocity, angle)
