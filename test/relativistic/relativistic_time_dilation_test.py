from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic import relativistic_time_dilation

# Using calculations from the paper: https://www.omnicalculator.com/physics/time-dilation

Args = namedtuple("Args", ["own_time", "velocity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    own_time = Quantity(10 * units.second)
    velocity = Quantity(200_000_000 * (units.meter / units.second))
    return Args(own_time=own_time, velocity=velocity)


def test_basic_time(test_args: Args) -> None:
    result = relativistic_time_dilation.calculate_relativistic_time(test_args.own_time,
        test_args.velocity)
    assert_equal(result, 13.423 * units.second)


def test_bad_time(test_args: Args) -> None:
    bt = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_time_dilation.calculate_relativistic_time(bt, test_args.velocity)
    with raises(TypeError):
        relativistic_time_dilation.calculate_relativistic_time(1, test_args.velocity)


def test_bad_velocity(test_args: Args) -> None:
    bv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_time_dilation.calculate_relativistic_time(test_args.own_time, bv)
    with raises(TypeError):
        relativistic_time_dilation.calculate_relativistic_time(test_args.own_time, 100)
