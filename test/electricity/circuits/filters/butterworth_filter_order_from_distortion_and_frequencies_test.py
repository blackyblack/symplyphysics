from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits.filters import butterworth_filter_order_from_distortion_and_frequencies as order_law

# Description
## The cutoff frequency is 1 hertz. The frequency of the band-stop is 5 Hertz. The distortion in the bandwidth is 0.5, the band-stop distortion
## is 31 and the Butterworth filter function will be used as an approximating function. Then the order of the Butterworth filter will be 3.
## https://ru.dsplib.org/content/filter_butter_ap/filter_butter_ap.html#r3

Args = namedtuple("Args",
    ["bandwidth_distortion", "band_stop_distortion", "cutoff_frequency", "band_stop_frequency"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    bandwidth_distortion = 0.5
    band_stop_distortion = 31
    cutoff_frequency = Quantity(1 * units.hertz)
    band_stop_frequency = Quantity(5 * units.hertz)
    return Args(bandwidth_distortion=bandwidth_distortion,
        band_stop_distortion=band_stop_distortion,
        band_stop_frequency=band_stop_frequency,
        cutoff_frequency=cutoff_frequency)


def test_basic_butterworth_filter_order(test_args: Args) -> None:
    result = order_law.calculate_butterworth_filter_order(test_args.bandwidth_distortion,
        test_args.band_stop_distortion, test_args.band_stop_frequency, test_args.cutoff_frequency)
    assert_equal(result, 3)


def test_bad_distortions(test_args: Args) -> None:
    bad_distortion = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_butterworth_filter_order(bad_distortion, test_args.band_stop_distortion,
            test_args.band_stop_frequency, test_args.cutoff_frequency)
    with raises(errors.UnitsError):
        order_law.calculate_butterworth_filter_order(test_args.bandwidth_distortion, bad_distortion,
            test_args.band_stop_frequency, test_args.cutoff_frequency)


def test_bad_frequencies(test_args: Args) -> None:
    bad_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_butterworth_filter_order(test_args.bandwidth_distortion,
            test_args.band_stop_distortion, bad_frequency, test_args.cutoff_frequency)
    with raises(TypeError):
        order_law.calculate_butterworth_filter_order(test_args.bandwidth_distortion,
            test_args.band_stop_distortion, 100, test_args.cutoff_frequency)
    with raises(errors.UnitsError):
        order_law.calculate_butterworth_filter_order(test_args.bandwidth_distortion,
            test_args.band_stop_distortion, test_args.band_stop_frequency, bad_frequency)
    with raises(TypeError):
        order_law.calculate_butterworth_filter_order(test_args.bandwidth_distortion,
            test_args.band_stop_distortion, test_args.band_stop_frequency, 100)
