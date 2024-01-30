from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.relativistic import coordinate_conversion_velocity_constant as coordinate_law

# Description
## Let the coordinate in the first frame of reference be 2 meter, the speed of the second system relative to
## the first is 3 meter per second, and the time_first_frame is 3 second.
## Then the coordinate in the second frame of reference will be -7 meter.
## https://physics.icalculator.com/lorentz-transformation-of-coordinates-calculator.html

@fixture(name="test_args")
def test_args_fixture():
    coordinate_first_frame = Quantity(2 * units.meter)
    velocity = Quantity(3 * (units.meter / units.second))
    time_first_frame = Quantity(3 * units.second)

    Args = namedtuple("Args", ["coordinate_first_frame", "velocity", "time_first_frame"])
    return Args(coordinate_first_frame=coordinate_first_frame, velocity=velocity, time_first_frame=time_first_frame)


def test_basic_another_coordinate(test_args):
    result = coordinate_law.calculate_coordinate_second_frame(test_args.coordinate_first_frame, test_args.velocity, test_args.time_first_frame)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result = convert_to(result, units.meter).evalf(5)
    assert result == approx(-7, 0.01)


def test_bad_coordinate(test_args):
    coordinate_first_frame = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        coordinate_law.calculate_coordinate_second_frame(coordinate_first_frame, test_args.velocity, test_args.time_first_frame)
    with raises(TypeError):
        coordinate_law.calculate_coordinate_second_frame(100, test_args.velocity, test_args.time_first_frame)


def test_bad_velocity(test_args):
    velocity = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        coordinate_law.calculate_coordinate_second_frame(test_args.coordinate_first_frame, velocity, test_args.time_first_frame)
    with raises(TypeError):
        coordinate_law.calculate_coordinate_second_frame(test_args.coordinate_first_frame, 100, test_args.time_first_frame)


def test_bad_time_first_frame(test_args):
    time_first_frame = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        coordinate_law.calculate_coordinate_second_frame(test_args.coordinate_first_frame, test_args.velocity, time_first_frame)
    with raises(TypeError):
        coordinate_law.calculate_coordinate_second_frame(test_args.coordinate_first_frame, test_args.velocity, 100)
