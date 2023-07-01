from collections import namedtuple
from pytest import approx, fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.relativistic.waves import longitudinal_frequency_shift_from_absolute_velocities as doppler_law

# Using calculations from the paper: http://www.mrelativity.net/TransDoppler/Relativistic%20Transverse%20Doppler%20Effect.pdf.


@fixture(name="test_args")
def test_args_fixture():
    # speed of light * 0.9
    object_velocity = Quantity(269813212.2 * units.meter / units.second)
    # emitted wavelength is 5.5 * 10^-7 meters (frequency is 5.45 * 10^14 Hz)
    emitted_frequency = Quantity(5.45077e14 * units.hertz)
    # observer is not moving so angle does not matter
    zero_velocity = Quantity(0 * units.meter / units.second)
    Args = namedtuple("Args", ["object_velocity", "emitted_frequency", "zero_velocity"])
    return Args(object_velocity=object_velocity,
        emitted_frequency=emitted_frequency,
        zero_velocity=zero_velocity)


def test_basic_frequency(test_args):
    result = doppler_law.calculate_observed_frequency(test_args.emitted_frequency, speed_of_light,
        test_args.object_velocity, test_args.zero_velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    result_freq = int(convert_to(result, units.hertz).subs(units.hertz, 1).evalf(4))
    assert result_freq == approx(1.2507e14, 0.001)


def test_classical_moving_observer_frequency(test_args):
    # observer is immobile and emitter is moving
    horn_frequency = Quantity(2000 * units.hertz)
    object_velocity = Quantity(-9 * units.kilometer / units.hour)
    # speed of sound
    wave_velocity = Quantity(340 * units.meter / units.second)
    result = doppler_law.calculate_observed_frequency(horn_frequency, wave_velocity,
        object_velocity, test_args.zero_velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    result_freq = int(convert_to(result, units.hertz).subs(units.hertz, 1).evalf(4))
    assert result_freq == approx(1985, 0.1)

    # make observer moving and source idle to verify that effect is irrelative if we have moving
    # source of waves or moving observer
    result = doppler_law.calculate_observed_frequency(horn_frequency, wave_velocity,
        test_args.zero_velocity, object_velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    moving_observer_freq = int(convert_to(result, units.hertz).subs(units.hertz, 1).evalf(4))
    assert moving_observer_freq == result_freq


def test_moving_observer_frequency(test_args):
    # observer is immobile and emitter is moving
    # speed of light * 0.8
    object_velocity = Quantity(239833966 * units.meter / units.second)
    result = doppler_law.calculate_observed_frequency(test_args.emitted_frequency, speed_of_light,
        object_velocity, test_args.zero_velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    result_freq = int(convert_to(result, units.hertz).subs(units.hertz, 1).evalf(8))

    # make observer moving and source idle to verify that effect is irrelative if we have moving
    # source of waves or moving observer
    result = doppler_law.calculate_observed_frequency(test_args.emitted_frequency, speed_of_light,
        test_args.zero_velocity, object_velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    moving_observer_freq = int(convert_to(result, units.hertz).subs(units.hertz, 1).evalf(8))
    assert moving_observer_freq == result_freq


def test_bad_velocity(test_args):
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, speed_of_light, vb,
            test_args.zero_velocity)
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, speed_of_light, 100,
            test_args.zero_velocity)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, speed_of_light,
            test_args.object_velocity, vb)
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, speed_of_light,
            test_args.object_velocity, 100)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, vb,
            test_args.object_velocity, test_args.zero_velocity)
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(test_args.emitted_frequency, 100,
            test_args.object_velocity, test_args.zero_velocity)


def test_bad_frequency(test_args):
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        doppler_law.calculate_observed_frequency(fb, speed_of_light, test_args.object_velocity,
            test_args.zero_velocity)
    with raises(TypeError):
        doppler_law.calculate_observed_frequency(100, speed_of_light, test_args.object_velocity,
            test_args.zero_velocity)
