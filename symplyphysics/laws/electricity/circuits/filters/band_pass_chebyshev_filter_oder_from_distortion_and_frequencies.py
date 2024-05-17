from math import ceil
from sympy import (Eq, solve, acosh)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, dimensionless,
    convert_to_float)

# Description
## The approximation of the power transmission coefficient of a normalized high-pass filter is given by approximating functions of the order of n.
## The Chebyshev filter is described by the function 1 / cos(n * acos(1 / frequency)).
## Bandwidth distortion determines the maximum distortion in the bandwidth. In other words, to determine the level of ripples in the bandwidth.
## The band-stop distortion sets the required suppression level in the filter band-stop.
## You can only implement an integer-order filter. And the formula itself returns the minimum order of the filter, which will have the
## specified parameters (frequencies and distortions). Therefore, you need to round it up (use ceil()).
## In this case, a band-pass filter is considered. The bandpass filter only passes frequencies from a certain frequency range.
## http://www.dsplib.ru/content/filters/ch8/ch8.html#r4

## Law is: N = acosh(es / ep) / acosh((fz^2 - f0^2) / (delta_w * fz)), where
## N - Chebyshev filter order,
## es - band-stop distortion,
## ep - bandwidth distortion,
## fz - frequency of the band-stop (the frequency that corresponds to the maximum band-stop distortion),
## f0 - cutoff frequency (the frequency that corresponds to the maximum bandwidth distortion),
## delta_w - bandwidth.

band_pass_chebyshev_filter_order = Symbol("band_pass_chebyshev_filter_order", dimensionless)

bandwidth_distortion = Symbol("bandwidth_distortion", dimensionless)
band_stop_distortion = Symbol("band_stop_distortion", dimensionless)
cutoff_frequency = Symbol("cutoff_frequency", units.frequency)
band_stop_frequency = Symbol("band_stop_frequency", units.frequency)
bandwidth = Symbol("bandwidth", units.frequency)

law = Eq(
    band_pass_chebyshev_filter_order,
    acosh(bandwidth_distortion / band_stop_distortion) / acosh(
    (band_stop_frequency**2 - cutoff_frequency**2) / (bandwidth * band_stop_frequency)))


@validate_input(bandwidth_distortion_=bandwidth_distortion,
    band_stop_distortion_=band_stop_distortion,
    band_stop_frequency_=band_stop_frequency,
    cutoff_frequency_=cutoff_frequency,
    bandwidth_=bandwidth)
@validate_output(band_pass_chebyshev_filter_order)
def calculate_band_pass_chebyshev_filter_order(bandwidth_distortion_: float,
    band_stop_distortion_: float, band_stop_frequency_: Quantity, cutoff_frequency_: Quantity,
    bandwidth_: Quantity) -> int:
    result_expr = solve(law, band_pass_chebyshev_filter_order,
        dict=True)[0][band_pass_chebyshev_filter_order]
    result_expr = result_expr.subs({
        bandwidth_distortion: bandwidth_distortion_,
        band_stop_distortion: band_stop_distortion_,
        band_stop_frequency: band_stop_frequency_,
        cutoff_frequency: cutoff_frequency_,
        bandwidth: bandwidth_,
    })
    return ceil(convert_to_float(result_expr))
