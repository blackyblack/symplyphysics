from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.damped_oscillations import q_factor_via_bandwidth as q_factor_law

# Description
## The resonant frequency of the oscillator is 100 Hz and its bandwidth is 5 Hz. 
## The Q factor of the oscillator is 20.

Args = namedtuple("Args", "f df")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = Quantity(100 * units.hertz)
    df = Quantity(5 * units.hertz)
    return Args(f=f, df=df)


def test_law(test_args: Args) -> None:
    result = q_factor_law.calculate_q_factor(test_args.f, test_args.df)
    assert_equal(result, 20)


def test_bad_frequency(test_args: Args) -> None:
    fb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        q_factor_law.calculate_q_factor(fb, test_args.df)
    with raises(TypeError):
        q_factor_law.calculate_q_factor(100, test_args.df)
    with raises(errors.UnitsError):
        q_factor_law.calculate_q_factor(test_args.f, fb)
    with raises(TypeError):
        q_factor_law.calculate_q_factor(test_args.f, 100)
