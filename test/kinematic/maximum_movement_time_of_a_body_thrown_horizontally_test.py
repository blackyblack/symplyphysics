from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.kinematic import maximum_movement_time_of_a_body_thrown_horizontally as time_law

# Description
## Let the height be 10 meter, and the acceleration of gravity is 9.8 meter per second.
## Then the drop time will be 1.428 second.
## https://www.omnicalculator.com/physics/horizontal-projectile-motion


@fixture(name="test_args")
def test_args_fixture():
    height = Quantity(10 * units.meter)
    acceleration = Quantity(9.8 * (units.meter / units.second**2))

    Args = namedtuple("Args", ["height", "acceleration"])
    return Args(height=height, acceleration=acceleration)


def test_basic_movement_time(test_args):
    result = time_law.calculate_movement_time(test_args.height, test_args.acceleration)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.time)
    result = convert_to(result, units.second).evalf(5)
    assert result == approx(1.428, rel=0.01)


def test_bad_height(test_args):
    height = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        time_law.calculate_movement_time(height, test_args.acceleration)
    with raises(TypeError):
        time_law.calculate_movement_time(100, test_args.acceleration)


def test_bad_acceleration(test_args):
    acceleration = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        time_law.calculate_movement_time(test_args.height, acceleration)
    with raises(TypeError):
        time_law.calculate_movement_time(test_args.height, 100)
