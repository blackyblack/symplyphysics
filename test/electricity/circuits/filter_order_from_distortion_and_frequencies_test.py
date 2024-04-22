from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity
)
from symplyphysics.laws.electricity.circuits import filter_order_from_distortion_and_frequencies as order_law

# Description
## The cutoff frequency is 1 hertz. The frequency of the band-stop is 5 Hertz. The distortion in the bandwidth is 0.5, the band-stop distortion
## is 31 and the Butterworth filter function will be used as an approximating function. Then the order of the Butterworth filter will be 3.
## https://ru.dsplib.org/content/filter_butter_ap/filter_butter_ap.html#r3

Args = namedtuple("Args", [
    "filter_function", "band_stop_frequency", "bandwidth_distortion", "band_stop_distortion"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    band_stop_frequency = Quantity(5 * units.hertz)
    cutoff_frequency = Quantity(1 * units.hertz)
    filter_function = (band_stop_frequency / cutoff_frequency)**order_law.filter_order
    bandwidth_distortion = 0.5
    band_stop_distortion = 31

    return Args(filter_function=filter_function,
        band_stop_frequency=band_stop_frequency,
        bandwidth_distortion=bandwidth_distortion,
        band_stop_distortion=band_stop_distortion)


def test_basic_order(test_args: Args) -> None:
    result = order_law.calculate_order(test_args.filter_function, test_args.band_stop_frequency, test_args.bandwidth_distortion, test_args.band_stop_distortion)
    assert_equal(result, 3, tolerance=0.01)


def test_bad_band_stop_frequency(test_args: Args) -> None:
    band_stop_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_order(test_args.filter_function, band_stop_frequency, test_args.bandwidth_distortion, test_args.band_stop_distortion)
    with raises(TypeError):
        order_law.calculate_order(test_args.filter_function, 100, test_args.bandwidth_distortion, test_args.band_stop_distortion)


def test_bad_distortions(test_args: Args) -> None:
    bad_distortion = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_order(test_args.filter_function, test_args.band_stop_frequency, bad_distortion, test_args.band_stop_distortion)
    with raises(errors.UnitsError):
        order_law.calculate_order(test_args.filter_function, test_args.band_stop_frequency, test_args.bandwidth_distortion, bad_distortion)
