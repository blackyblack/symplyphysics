from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import centripetal_acceleration_is_squared_velocity_by_radius as centripetal_acceleration_def

# Description
## The object is rotating with linear velocity of 10 m/s tied to the center with 0/5 m thread, it's centripetal acceleration should be 200 m/s/s (according to calc.ru online calculator)


@fixture
def test_args():
    lin_velocity = Quantity(10 * units.meter / units.second)
    curve_radius = Quantity(0.5 * units.meter)
    Args = namedtuple("Args", ["lin_velocity", "curve_radius"])
    return Args(lin_velocity=lin_velocity, curve_radius=curve_radius)


def test_basic_acceleration(test_args):
    result = centripetal_acceleration_def.calculate_acceleration(
        test_args.lin_velocity, test_args.curve_radius)
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
                                                     units.acceleration)
    result_acceleration = convert_to(result,
                                     units.meter / units.second**2).subs({
                                         units.meter: 1,
                                         units.second: 1
                                     }).evalf(2)
    assert result_acceleration == approx(200.0, 0.01)


def test_bad_velocity(test_args):
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        centripetal_acceleration_def.calculate_acceleration(
            vb, test_args.curve_radius)
    with raises(TypeError):
        centripetal_acceleration_def.calculate_acceleration(
            100, test_args.curve_radius)


def test_bad_radius(test_args):
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        centripetal_acceleration_def.calculate_acceleration(
            test_args.lin_velocity, rb)
    with raises(TypeError):
        centripetal_acceleration_def.calculate_acceleration(
            test_args.lin_velocity, 100)
