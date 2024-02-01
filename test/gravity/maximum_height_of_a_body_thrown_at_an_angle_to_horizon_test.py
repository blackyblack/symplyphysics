from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.gravity import maximum_height_of_a_body_thrown_at_an_angle_to_horizon as height_law

# Description
## Let the initial velocity be 3 meter per second, the angle of the throw is 45 degree (pi / 4 radian),
## and the acceleration of gravity is 9.8 [meter / second^2]. Then the maximum lifting height of the body thrown at
## an angle to the horizon will be 0.23 meter.
## https://www.omnicalculator.com/physics/projectile-motion


@fixture(name="test_args")
def test_args_fixture():
    initial_velocity = Quantity(3 * (units.meter / units.second))
    angle = pi / 4

    Args = namedtuple("Args", ["initial_velocity", "angle"])
    return Args(initial_velocity=initial_velocity, angle=angle)


def test_basic_height(test_args):
    result = height_law.calculate_height(test_args.initial_velocity, test_args.angle)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result = convert_to(result, units.meter).evalf(5)
    assert result == approx(0.23, rel=0.01)


def test_bad_initial_velocity(test_args):
    initial_velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        height_law.calculate_height(initial_velocity, test_args.angle)
    with raises(TypeError):
        height_law.calculate_height(100, test_args.angle)


def test_bad_angle(test_args):
    angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        height_law.calculate_height(test_args.initial_velocity, angle)
