from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.kinematic import maximum_movement_time_of_a_body_thrown_at_an_angle_to_horizon as time_law

# Description
## Let the initial velocity be 3 meter per second, the angle of the throw is 45 degree (pi / 4 radian),
## and the acceleration of gravity is 9.8 [meter / second^2]. Then the travel time will be 0.43 second.
## https://www.indigomath.ru//raschety/yJTxXc.html

@fixture(name="test_args")
def test_args_fixture():
    initial_velocity = Quantity(3 * (units.meter / units.second))
    angle = pi / 4
    acceleration = Quantity(9.8 * (units.meter / units.second**2))

    Args = namedtuple("Args", ["initial_velocity", "angle", "acceleration"])
    return Args(initial_velocity=initial_velocity, angle=angle, acceleration=acceleration)


def test_basic_movement_time(test_args):
    result = time_law.calculate_movement_time(test_args.initial_velocity, test_args.angle,
        test_args.acceleration)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.time)
    result = convert_to(result, units.second).evalf(5)
    assert result == approx(0.43, rel=0.01)


def test_bad_initial_velocity(test_args):
    initial_velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_movement_time(initial_velocity, test_args.angle,
            test_args.acceleration)
    with raises(TypeError):
        time_law.calculate_movement_time(100, test_args.angle,
            test_args.acceleration)


def test_bad_angle(test_args):
    angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_movement_time(test_args.initial_velocity, angle,
            test_args.acceleration)
    with raises(AttributeError):
        time_law.calculate_movement_time(test_args.initial_velocity, True, test_args.acceleration)


def test_bad_acceleration(test_args):
    acceleration = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        time_law.calculate_movement_time(test_args.initial_velocity, test_args.angle,
            acceleration)
    with raises(TypeError):
        time_law.calculate_movement_time(test_args.initial_velocity, test_args.angle, 100)
