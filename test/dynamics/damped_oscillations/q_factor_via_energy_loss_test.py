from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.damped_oscillations import q_factor_via_energy_loss as q_factor_law

# Description
## An oscillator with resonant angular frequency w = 100 rad/s has 5 J of energy stored in it.
## The power loss of the oscillator is 200 W. The Q factor of the oscillator amounts to 2.5.

Args = namedtuple("Args", "w e p")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w = Quantity(100 * units.radian / units.second)
    e = Quantity(5 * units.joule)
    p = Quantity(200 * units.watt)
    return Args(w=w, e=e, p=p)


def test_law(test_args: Args) -> None:
    result = q_factor_law.calculate_q_factor(test_args.w, test_args.e, test_args.p)
    assert_equal(result, 2.5)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        q_factor_law.calculate_q_factor(wb, test_args.e, test_args.p)
    with raises(TypeError):
        q_factor_law.calculate_q_factor(100, test_args.e, test_args.p)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        q_factor_law.calculate_q_factor(test_args.w, eb, test_args.p)
    with raises(TypeError):
        q_factor_law.calculate_q_factor(test_args.w, 100, test_args.p)


def test_bad_power(test_args: Args) -> None:
    pb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        q_factor_law.calculate_q_factor(test_args.w, test_args.e, pb)
    with raises(TypeError):
        q_factor_law.calculate_q_factor(test_args.w, test_args.e, 100)
