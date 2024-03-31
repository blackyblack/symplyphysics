from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic.waves import longitudinal_frequency_shift_from_velocity as doppler_law

# Description. With help of online calculator at https://physics.icalculator.com/the-doppler-effect-in-light-waves-calculator.html:
## For red wave with frequency 384THz (3.84e14 Hz) emitted by object moving away from observer with velocity 0.1c = 29979245.8 m/s
## observed frequency should be 347341068940800

Args = namedtuple("Args", ["object_velocity", "emitted_frequency"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    object_velocity = Quantity(29979245.8 * units.meter / units.second)
    emitted_frequency = Quantity(3.84e14 * units.hertz)
    return Args(object_velocity=object_velocity, emitted_frequency=emitted_frequency)


def test_basic_frequency(test_args: Args) -> None:
    result = doppler_law.calculate_observed_frequency(test_args.emitted_frequency,
        test_args.object_velocity)
    assert_equal(result, 3.473e14 * units.hertz)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, vb)
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, 100)


def test_bad_frequency(test_args: Args) -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(fb, test_args.object_velocity)
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(100, test_args.object_velocity)
