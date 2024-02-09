from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import period_from_angular_frequency as period_def

Args = namedtuple("Args", ["w"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w = Quantity(6.28 * units.radian / units.second)
    return Args(w=w)


def test_basic_period(test_args: Args) -> None:
    result = period_def.calculate_period(test_args.w)
    assert_equal(result, 1 * units.second)


def test_bad_frequency() -> None:
    wb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        period_def.calculate_period(wb)
    with raises(TypeError):
        period_def.calculate_period(100)
