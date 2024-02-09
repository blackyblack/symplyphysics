from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import period_of_ideal_pendulum_from_length as pendulum_period

Args = namedtuple("Args", ["L"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    L = Quantity(1 * units.meter)
    return Args(L=L)


def test_basic_period(test_args: Args) -> None:
    result = pendulum_period.calculate_period(test_args.L)
    # For a pendulum of 1 meter length, period should be 2 seconds
    assert_equal(result, 2.006 * units.second)


def test_bad_length() -> None:
    Lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pendulum_period.calculate_period(Lb)
    with raises(TypeError):
        pendulum_period.calculate_period(100)
