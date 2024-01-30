from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.relativistic import coordinate_conversion as coordinate_law

# Description
## Let the coordinate in the first frame of reference be 2 meter, the speed of the second system relative to
## the first is 3 meter per second, and the time is 3 second.
## Then the coordinate in the second frame of reference will be -7 meter.
## https://physics.icalculator.com/lorentz-transformation-of-coordinates-calculator.html

@fixture(name="test_args")
def test_args_fixture():
    coordinate = Quantity(2 * units.meter)
    velocity = Quantity(3 * (units.meter / units.second))
    time = Quantity(3 * units.second)

    Args = namedtuple("Args", ["coordinate", "velocity", "time"])
    return Args(coordinate=coordinate, velocity=velocity, time=time)


def test_basic_another_coordinate(test_args):
    result = coordinate_law.calculate_another_coordinate(test_args.coordinate, test_args.velocity, test_args.time)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_voltage = convert_to(result, units.meter).evalf(5)
    assert result_voltage == approx(-7, 0.01)


def test_bad_coordinate(test_args):
    coordinate = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        coordinate_law.calculate_another_coordinate(coordinate, test_args.velocity, test_args.time)
    with raises(TypeError):
        coordinate_law.calculate_another_coordinate(100, test_args.velocity, test_args.time)


def test_bad_velocity(test_args):
    velocity = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        coordinate_law.calculate_another_coordinate(test_args.coordinate, velocity, test_args.time)
    with raises(TypeError):
        coordinate_law.calculate_another_coordinate(test_args.coordinate, 100, test_args.time)


def test_bad_time(test_args):
    time = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        coordinate_law.calculate_another_coordinate(test_args.coordinate, test_args.velocity, time)
    with raises(TypeError):
        coordinate_law.calculate_another_coordinate(test_args.coordinate, test_args.velocity, 100)
