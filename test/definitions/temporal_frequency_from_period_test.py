from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import temporal_frequency_from_period as frequency_def

Args = namedtuple("Args", ["T"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    T = Quantity(2.8 * units.second)
    return Args(T=T)


def test_basic_frequency(test_args: Args) -> None:
    result = frequency_def.calculate_frequency(test_args.T)
    assert_equal(result, 0.3571 * units.hertz)


def test_bad_period() -> None:
    Tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(Tb)
    with raises(TypeError):
        frequency_def.calculate_frequency(100)
