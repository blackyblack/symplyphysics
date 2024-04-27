from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import band_pass_chebyshev_filter_oder_from_distortion_and_frequencies as order_law

# Description
## The cutoff frequency is 1 hertz. The upper stop frequency is 18 Hertz. The distortion in the bandwidth is 316, the band-stop distortion
## is 0.76 and the Chebyshev filter function will be used as an approximating function. The bandwidth is 10 herz. Then the order of the Chebyshev filter will be 6.

Args = namedtuple("Args", ["bandwidth_distortion", "band_stop_distortion", "cutoff_frequency", "band_stop_frequency", "bandwidth"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    bandwidth_distortion = 316
    band_stop_distortion = 0.76
    cutoff_frequency = Quantity(1 * units.hertz)
    band_stop_frequency = Quantity(18 * units.hertz)
    bandwidth = Quantity(10 * units.hertz)
    return Args(bandwidth_distortion=bandwidth_distortion,
        band_stop_distortion=band_stop_distortion,
        band_stop_frequency=band_stop_frequency,
        cutoff_frequency=cutoff_frequency,
        bandwidth=bandwidth)


def test_basic_band_pass_chebyshev_filter_order(test_args: Args) -> None:
    result = order_law.calculate_band_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion,
        test_args.band_stop_frequency, test_args.cutoff_frequency, test_args.bandwidth)
    assert_equal(result, 6)


def test_bad_distortions(test_args: Args) -> None:
    bad_distortion = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_band_pass_chebyshev_filter_order(bad_distortion, test_args.band_stop_distortion, test_args.band_stop_frequency,
            test_args.cutoff_frequency, test_args.bandwidth)
    with raises(errors.UnitsError):
        order_law.calculate_band_pass_chebyshev_filter_order(test_args.bandwidth_distortion, bad_distortion, test_args.band_stop_frequency,
            test_args.cutoff_frequency, test_args.bandwidth)


def test_bad_frequencies(test_args: Args) -> None:
    bad_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_band_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, bad_frequency,
            test_args.cutoff_frequency, test_args.bandwidth)
    with raises(TypeError):
        order_law.calculate_band_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, 100,
            test_args.cutoff_frequency, test_args.bandwidth)
    with raises(errors.UnitsError):
        order_law.calculate_band_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, test_args.band_stop_frequency,
            bad_frequency, test_args.bandwidth)
    with raises(TypeError):
        order_law.calculate_band_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, test_args.band_stop_frequency,
            100, test_args.bandwidth)


def test_bad_bandwidth(test_args: Args) -> None:
    bad_bandwidth = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        order_law.calculate_band_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, test_args.band_stop_frequency,
            test_args.cutoff_frequency, bad_bandwidth)
    with raises(TypeError):
        order_law.calculate_band_pass_chebyshev_filter_order(test_args.bandwidth_distortion, test_args.band_stop_distortion, test_args.band_stop_frequency,
            test_args.cutoff_frequency, 100)
