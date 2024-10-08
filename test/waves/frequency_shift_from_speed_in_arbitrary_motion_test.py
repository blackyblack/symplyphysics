from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    angle_type,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import frequency_shift_from_speed_in_arbitrary_motion as doppler_law

## Man stands near railroad and hears horn of incoming train. Angle between line of sight to the train is 30 degrees.
## Train comes with velocity of 9km/h and horns with 2000Hz frequency. What frequency hears the man?
## We have online calc for Doppler effect here: https://planetcalc.ru/2351/. With our parameters we should obtain 2013Hz observed frequency.
## Another situation is when man rides a bike towards horning standing train with the same velocity. Observed frequency should be same.

Args = namedtuple("Args", [
    "sound_velocity", "train_speed", "bike_speed", "zero_velocity", "horn_frequency",
    "source_angle", "zero_angle"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    sound_velocity = Quantity(340 * units.meter / units.second)
    train_speed = Quantity(9 * units.kilometer / units.hour)
    bike_speed = Quantity(9 * units.kilometer / units.hour)
    zero_velocity = Quantity(0)
    horn_frequency = Quantity(2000 * units.hertz)
    source_angle = Quantity(30 * units.degree, dimension=angle_type)
    zero_angle = Quantity(0 * units.radian, dimension=angle_type)
    return Args(sound_velocity=sound_velocity,
        zero_velocity=zero_velocity,
        train_speed=train_speed,
        bike_speed=bike_speed,
        horn_frequency=horn_frequency,
        source_angle=source_angle,
        zero_angle=zero_angle)


def test_basic_frequency(test_args: Args) -> None:
    result_1 = doppler_law.calculate_observed_frequency(test_args.horn_frequency,
        test_args.sound_velocity, (test_args.train_speed, test_args.source_angle),
        (test_args.zero_velocity, test_args.zero_angle))
    # Observer is moving towards the source of sound, hence the angle between signal vector and his
    # direction is over 90 degrees
    observer_angle = Quantity((180 - 30) * units.degree, dimension=angle_type)
    result_2 = doppler_law.calculate_observed_frequency(test_args.horn_frequency,
        test_args.sound_velocity, (test_args.zero_velocity, test_args.zero_angle),
        (test_args.bike_speed, observer_angle))
    assert_equal(result_1, 2013 * units.hertz)
    # Doppler effect is irrelative at relatively low velocities
    assert_equal(result_2, 2013 * units.hertz)


# Classical Doppler effect does not have frequency shift when moving at 90 degrees
def test_transverse_frequency(test_args: Args) -> None:
    source_angle = Quantity(90 * units.degree, dimension=angle_type)
    result = doppler_law.calculate_observed_frequency(test_args.horn_frequency,
        test_args.sound_velocity, (test_args.train_speed, source_angle),
        (test_args.zero_velocity, test_args.zero_angle))
    assert_equal(result, test_args.horn_frequency)

    result = doppler_law.calculate_observed_frequency(test_args.horn_frequency,
        test_args.sound_velocity, (test_args.zero_velocity, test_args.zero_angle),
        (test_args.train_speed, source_angle))
    assert_equal(result, test_args.horn_frequency)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.horn_frequency, vb,
            (test_args.train_speed, test_args.source_angle),
            (test_args.bike_speed, test_args.zero_angle))
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(test_args.horn_frequency, 100,
            (test_args.train_speed, test_args.source_angle),
            (test_args.bike_speed, test_args.zero_angle))
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.horn_frequency, test_args.sound_velocity,
            (vb, test_args.source_angle), (test_args.bike_speed, test_args.zero_angle))
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(test_args.horn_frequency, test_args.sound_velocity,
            (100, test_args.source_angle), (test_args.bike_speed, test_args.zero_angle))
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.horn_frequency, test_args.sound_velocity,
            (test_args.train_speed, test_args.source_angle), (vb, test_args.zero_angle))
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(test_args.horn_frequency, test_args.sound_velocity,
            (test_args.train_speed, test_args.source_angle), (100, test_args.zero_angle))


def test_bad_frequency(test_args: Args) -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(fb, test_args.sound_velocity,
            (test_args.train_speed, test_args.source_angle),
            (test_args.bike_speed, test_args.zero_angle))
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(100, test_args.sound_velocity,
            (test_args.train_speed, test_args.source_angle),
            (test_args.bike_speed, test_args.zero_angle))


def test_bad_angle(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.horn_frequency, test_args.sound_velocity,
            (test_args.train_speed, ab), (test_args.bike_speed, test_args.zero_angle))
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.horn_frequency, test_args.sound_velocity,
            (test_args.train_speed, test_args.source_angle), (test_args.bike_speed, ab))
