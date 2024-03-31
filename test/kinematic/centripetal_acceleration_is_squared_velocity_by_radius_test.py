from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import centripetal_acceleration_is_squared_velocity_by_radius as centripetal_acceleration_def

# Description
## The object is rotating with linear velocity of 10 m/s tied to the center with 0/5 m thread, it's centripetal acceleration should be 200 m/s/s (according to calc.ru online calculator)

Args = namedtuple("Args", ["lin_velocity", "curve_radius"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    lin_velocity = Quantity(10 * units.meter / units.second)
    curve_radius = Quantity(0.5 * units.meter)
    return Args(lin_velocity=lin_velocity, curve_radius=curve_radius)


def test_basic_acceleration(test_args: Args) -> None:
    result = centripetal_acceleration_def.calculate_acceleration(test_args.lin_velocity,
        test_args.curve_radius)
    assert_equal(result, 200 * units.meter / units.second**2)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        centripetal_acceleration_def.calculate_acceleration(vb, test_args.curve_radius)
    with raises(TypeError):
        centripetal_acceleration_def.calculate_acceleration(100, test_args.curve_radius)


def test_bad_radius(test_args: Args) -> None:
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        centripetal_acceleration_def.calculate_acceleration(test_args.lin_velocity, rb)
    with raises(TypeError):
        centripetal_acceleration_def.calculate_acceleration(test_args.lin_velocity, 100)
