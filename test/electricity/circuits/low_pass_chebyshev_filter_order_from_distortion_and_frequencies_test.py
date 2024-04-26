from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import low_pass_chebyshev_filter_order_from_distortion_and_frequencies as order_law

# Description
## The cutoff frequency is 0.48 hertz. The frequency of the band-stop is 1 Hertz. The distortion in the bandwidth is 0.76, the band-stop distortion
## is 316 and the Chebyshev filter function will be used as an approximating function. Then the order of the Chebyshev filter will be 5.
## https://ru.dsplib.org/content/filter_cheby2_ap/filter_cheby2_ap.html#r3

Args = namedtuple("Args", ["bandwidth_distortion", "band_stop_distortion", "cutoff_frequency", "band_stop_frequency"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    bandwidth_distortion = 0.76
    band_stop_distortion = 316
    cutoff_frequency = Quantity(0.48 * units.hertz)
    band_stop_frequency = Quantity(1 * units.hertz)
    return Args(bandwidth_distortion=bandwidth_distortion,
        band_stop_distortion=band_stop_distortion,
        band_stop_frequency=band_stop_frequency,
        cutoff_frequency=cutoff_frequency)


def test_basic_low_pass_chebyshev_filter_order(test_args: Args) -> None:
    result = order_law.calculate_low_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion,
        test_args.band_stop_frequency, test_args.cutoff_frequency)
    assert_equal(result, 5)


def test_bad_distortions(test_args: Args) -> None:
    bad_distortion = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_low_pass_chebyshev_filter_order(bad_distortion, test_args.band_stop_distortion, test_args.band_stop_frequency,
            test_args.cutoff_frequency)
    with raises(errors.UnitsError):
        order_law.calculate_low_pass_chebyshev_filter_order(test_args.bandwidth_distortion, bad_distortion, test_args.band_stop_frequency,
            test_args.cutoff_frequency)


def test_bad_frequencies(test_args: Args) -> None:
    bad_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_low_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, bad_frequency,
            test_args.cutoff_frequency)
    with raises(TypeError):
        order_law.calculate_low_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, 100,
            test_args.cutoff_frequency)
    with raises(errors.UnitsError):
        order_law.calculate_low_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, test_args.band_stop_frequency,
            bad_frequency)
    with raises(TypeError):
        order_law.calculate_low_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, test_args.band_stop_frequency,
            100)
