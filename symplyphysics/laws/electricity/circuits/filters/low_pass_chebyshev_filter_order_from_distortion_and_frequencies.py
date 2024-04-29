from math import ceil
from sympy import (Eq, solve, acos)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## The approximation of the power transmission coefficient of a normalized low-pass filter is given by approximating functions of the order of n.
## The Chebyshev filter is described by the function 1 / cos(n * acos(1 / frequency)).
## Bandwidth distortion determines the maximum distortion in the bandwidth. In other words, to determine the level of ripples in the bandwidth.
## The band-stop distortion sets the required suppression level in the filter band-stop.
## In this case, a low-pass filter is considered. The low-pass filter passes all frequencies from 0 to the set frequency.

## Law is: N = acos(es / ep) / acos(fs / fp), where
## N - low-pass Chebyshev filter order,
## es - band-stop distortion,
## ep - bandwidth distortion,
## fs - frequency of the band-stop (the frequency that corresponds to the maximum band-stop distortion),
## fp - cutoff frequency (the frequency that corresponds to the maximum bandwidth distortion).

low_pass_chebyshev_filter_order = Symbol("low_pass_chebyshev_filter_order", dimensionless)

bandwidth_distortion = Symbol("bandwidth_distortion", dimensionless)
band_stop_distortion = Symbol("band_stop_distortion", dimensionless)
cutoff_frequency = Symbol("cutoff_frequency", units.frequency)
band_stop_frequency = Symbol("band_stop_frequency", units.frequency)

law = Eq(low_pass_chebyshev_filter_order, acos(band_stop_distortion / bandwidth_distortion) / acos(band_stop_frequency / cutoff_frequency))


@validate_input(bandwidth_distortion_=bandwidth_distortion,
    band_stop_distortion_=band_stop_distortion,
    band_stop_frequency_=band_stop_frequency,
    cutoff_frequency_=cutoff_frequency)
@validate_output(low_pass_chebyshev_filter_order)
def calculate_low_pass_chebyshev_filter_order(bandwidth_distortion_: float, band_stop_distortion_: float, band_stop_frequency_: Quantity,
    cutoff_frequency_: Quantity) -> int:
    result_expr = solve(law, low_pass_chebyshev_filter_order, dict=True)[0][low_pass_chebyshev_filter_order]
    result_expr = result_expr.subs({
        bandwidth_distortion: bandwidth_distortion_,
        band_stop_distortion: band_stop_distortion_,
        band_stop_frequency: band_stop_frequency_,
        cutoff_frequency: cutoff_frequency_,
    })
    return ceil(convert_to_float(result_expr))
