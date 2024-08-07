from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import temporal_frequency_is_number_of_events_per_unit_time as frequency_def

# Description
## Pendulum is making 25 complete oscillations in 60 seconds. It should have frequency of 0.4167 Hz.

Args = namedtuple("Args", ["N", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    oscillations = 25
    t = Quantity(60 * units.second)
    return Args(N=oscillations, t=t)


def test_basic_frequency(test_args: Args) -> None:
    result = frequency_def.calculate_frequency(test_args.N, test_args.t)
    assert_equal(result, 0.4167 * units.hertz)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(test_args.N, tb)
    with raises(TypeError):
        frequency_def.calculate_frequency(test_args.N, 100)
