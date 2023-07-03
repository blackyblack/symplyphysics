from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import temporal_frequency_is_events_per_time as frequency_def

# Description
## Pendulum is making 25 complete oscillations in 60 seconds. It should have frequency of 0.416 Hz.


@fixture(name="test_args")
def test_args_fixture():
    oscillations = 25
    t = Quantity(60 * units.second)
    Args = namedtuple("Args", ["N", "t"])
    return Args(N=oscillations, t=t)


def test_basic_frequency(test_args):
    result = frequency_def.calculate_frequency(test_args.N, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    result_frequency = convert_to(result, frequency_def.definition_units_SI).evalf(2)
    assert result_frequency == approx(0.416, 0.01)


def test_bad_time(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(test_args.N, tb)
    with raises(TypeError):
        frequency_def.calculate_frequency(test_args.N, 100)
