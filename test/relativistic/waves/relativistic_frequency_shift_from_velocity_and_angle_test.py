from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    angle_type,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic.waves import frequency_shift_from_velocity_and_angle as doppler_law

# Using calculations from the paper: http://www.mrelativity.net/TransDoppler/Relativistic%20Transverse%20Doppler%20Effect.pdf.

Args = namedtuple("Args", ["object_velocity", "emitted_frequency", "source_angle"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # speed of light * 0.9
    object_velocity = Quantity(269813212.2 * units.meter / units.second)
    # emitted wavelength is 5.5 * 10^-7 meters (frequency is 5.45 * 10^14 Hz)
    emitted_frequency = Quantity(5.45077e14 * units.hertz)
    # angle of approach is 30 degrees (angle of recession is 150 degrees)
    source_angle = Quantity(30 * units.degree, dimension=angle_type)
    return Args(object_velocity=object_velocity,
        emitted_frequency=emitted_frequency,
        source_angle=source_angle)


def test_basic_frequency(test_args: Args) -> None:
    result = doppler_law.calculate_observed_frequency(test_args.emitted_frequency,
        test_args.object_velocity, test_args.source_angle)
    assert_equal(result, 1.0772e15 * units.hertz)


# Relativistic Doppler effect has visible frequency shift when moving at 90 degrees
def test_transverse_frequency(test_args: Args) -> None:
    source_angle = Quantity(90 * units.degree, dimension=angle_type)
    result = doppler_law.calculate_observed_frequency(test_args.emitted_frequency,
        test_args.object_velocity, source_angle)
    assert_equal(result, 2.3755e14 * units.hertz)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, vb,
            test_args.source_angle)
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, 100,
            test_args.source_angle)


def test_bad_frequency(test_args: Args) -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(fb, test_args.object_velocity,
            test_args.source_angle)
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(100, test_args.object_velocity,
            test_args.source_angle)


def test_bad_angle(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency,
            test_args.object_velocity, ab)
