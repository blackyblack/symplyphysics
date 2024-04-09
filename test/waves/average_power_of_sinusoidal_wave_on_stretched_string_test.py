from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.waves import average_power_of_sinusoidal_wave_on_stretched_string as power_law

# Description
## A wave is traveling with a phase speed v = 10 m/s along a stretched string of linear density mu = 10 g/m.
## The amplitude of the wave is 10 cm and the angular frequency of the wave is 8 rad/s. The average power
## of the wave is 32 mW.

Args = namedtuple("Args", "mu v w a")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mu = Quantity(10 * units.gram / units.meter)
    v = Quantity(10 * units.meter / units.second)
    w = Quantity(8 * units.radian / units.second)
    a = Quantity(10 * units.centimeter)
    return Args(mu=mu, v=v, w=w, a=a)


def test_law(test_args: Args) -> None:
    result = power_law.calculate_average_power(test_args.mu, test_args.v, test_args.w, test_args.a)
    assert_equal(result, 32 * prefixes.milli * units.watt)


def test_bad_linear_density(test_args: Args) -> None:
    mub = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_average_power(mub, test_args.v, test_args.w, test_args.a)
    with raises(TypeError):
        power_law.calculate_average_power(100, test_args.v, test_args.w, test_args.a)


def test_bad_phase_velocity(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_average_power(test_args.mu, vb, test_args.w, test_args.a)
    with raises(TypeError):
        power_law.calculate_average_power(test_args.mu, 100, test_args.w, test_args.a)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_average_power(test_args.mu, test_args.v, wb, test_args.a)
    with raises(TypeError):
        power_law.calculate_average_power(test_args.mu, test_args.v, 100, test_args.a)


def test_bad_amplitude(test_args: Args) -> None:
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_average_power(test_args.mu, test_args.v, test_args.w, ab)
    with raises(TypeError):
        power_law.calculate_average_power(test_args.mu, test_args.v, test_args.w, 100)
