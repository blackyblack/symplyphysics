from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import temporal_frequency_from_period as frequency_def


@fixture
def test_args():
    T = Quantity(2.8 * units.second)
    Args = namedtuple("Args", ["T"])
    return Args(T=T)


def test_basic_frequency(test_args):
    result = frequency_def.calculate_frequency(test_args.T)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    result_period = convert_to(result, units.hertz).subs(units.hertz, 1).evalf(2)
    assert result_period == approx(0.36, 0.01)


def test_bad_period():
    Tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(Tb)
    with raises(TypeError):
        frequency_def.calculate_frequency(100)
