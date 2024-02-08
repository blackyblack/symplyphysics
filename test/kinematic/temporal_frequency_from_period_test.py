from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import temporal_frequency_from_period as frequency_def


@fixture(name="test_args")
def test_args_fixture():
    T = Quantity(2.8 * units.second)
    Args = namedtuple("Args", ["T"])
    return Args(T=T)


def test_basic_frequency(test_args):
    result = frequency_def.calculate_frequency(test_args.T)
    assert_equal(result, 0.3571 * units.hertz)


def test_bad_period():
    Tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(Tb)
    with raises(TypeError):
        frequency_def.calculate_frequency(100)
