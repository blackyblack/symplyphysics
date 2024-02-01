from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.gravity import maximum_movement_time_of_a_body_thrown_horizontally as time_law

# Description
## Let the height be 10 meter, and the acceleration of gravity is 9.8 meter per second.
## Then the drop time will be 1.428 second.
## https://www.omnicalculator.com/physics/horizontal-projectile-motion


@fixture(name="test_args")
def test_args_fixture():
    height = Quantity(10 * units.meter)

    Args = namedtuple("Args", ["height"])
    return Args(height=height)


def test_basic_movement_time(test_args):
    result = time_law.calculate_movement_time(test_args.height)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.time)
    result = convert_to(result, units.second).evalf(5)
    assert result == approx(1.428, rel=0.01)


def test_bad_height():
    height = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        time_law.calculate_movement_time(height)
    with raises(TypeError):
        time_law.calculate_movement_time(100)
