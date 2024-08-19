from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    angle_type,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematics import average_angular_speed_is_angular_distance_over_time as frequency_def

# Description
## Circle is rotating to 9.42 radians in 5.4 seconds. It should have angular frequency of 1.744 rad/s.

Args = namedtuple("Args", ["N", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rotation = Quantity(9.42 * units.radian, dimension=angle_type)
    t = Quantity(5.4 * units.second)
    return Args(N=rotation, t=t)


def test_basic_frequency(test_args: Args) -> None:
    result = frequency_def.calculate_frequency(test_args.N, test_args.t)
    assert_equal(result, 1.744 * units.radian / units.second)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(test_args.N, tb)
    with raises(TypeError):
        frequency_def.calculate_frequency(test_args.N, 100)


def test_bad_radians(test_args: Args) -> None:
    Nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(Nb, test_args.t)
