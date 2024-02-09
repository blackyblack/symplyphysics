from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration

Args = namedtuple("Args", ["v0", "v1", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v0 = Quantity(1 * units.meter / units.second)
    v1 = Quantity(20 * units.meter / units.second)
    t = Quantity(5 * units.second)
    return Args(v0=v0, v1=v1, t=t)


def test_basic_acceleration(test_args: Args) -> None:
    result = acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, test_args.t)
    assert_equal(result, 3.8 * units.meter / units.second**2)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(vb, test_args.v1, test_args.t)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(test_args.v0, vb, test_args.t)
    with raises(TypeError):
        acceleration.calculate_linear_acceleration(100, test_args.v1, test_args.t)
    with raises(TypeError):
        acceleration.calculate_linear_acceleration(test_args.v0, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, tb)
    with raises(TypeError):
        acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, 100)
