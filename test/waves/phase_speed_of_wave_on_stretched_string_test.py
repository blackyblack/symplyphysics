from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import phase_speed_of_wave_on_stretched_string as speed_law

# Description
## A string of linear density mu = 10 g/m is stretched so that its tension tau is 1 N.
## The speed of the wave traveling on the string is 10 m/s.

Args = namedtuple("Args", "tau mu")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    tau = Quantity(1.0 * units.newton)
    mu = Quantity(10.0 * units.gram / units.meter)
    return Args(tau=tau, mu=mu)


def test_law(test_args: Args) -> None:
    result = speed_law.calculate_wave_speed(test_args.tau, test_args.mu)
    assert_equal(result, 10.0 * units.meter / units.second)


def test_bad_tension(test_args: Args) -> None:
    tau_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_wave_speed(tau_bad, test_args.mu)
    with raises(TypeError):
        speed_law.calculate_wave_speed(100, test_args.mu)


def test_bad_density(test_args: Args) -> None:
    mu_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_wave_speed(test_args.tau, mu_bad)
    with raises(TypeError):
        speed_law.calculate_wave_speed(test_args.tau, 100)
