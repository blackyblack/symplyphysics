from math import ceil
from sympy import (Eq, solve, log)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## The approximation of the power transmission coefficient of a normalized low-pass filter is given by approximating functions of the order of n.
## The Butterworth filter is described by the function frequency^n.
## Bandwidth distortion determines the maximum distortion in the bandwidth. In other words, to determine the level of ripples in the bandwidth.
## The band-stop distortion sets the required suppression level in the filter band-stop.

## Law is: N = ln(es / ep) / ln(fs / fp), where
## N - Butterworth filter order,
## es - band-stop distortion,
## ep - bandwidth distortion,
## fs - frequency of the band-stop (the frequency that corresponds to the maximum band-stop distortion),
## fp - cutoff frequency (the frequency that corresponds to the maximum bandwidth distortion).

butterworth_filter_order = Symbol("butterworth_filter_order", dimensionless)

bandwidth_distortion = Symbol("bandwidth_distortion", dimensionless)
band_stop_distortion = Symbol("band_stop_distortion", dimensionless)
cutoff_frequency = Symbol("cutoff_frequency", units.frequency)
band_stop_frequency = Symbol("band_stop_frequency", units.frequency)

law = Eq(butterworth_filter_order, log(band_stop_distortion / bandwidth_distortion) / log(band_stop_frequency / cutoff_frequency))


def print_law() -> str:
    return print_expression(law)


@validate_input(bandwidth_distortion_=bandwidth_distortion,
    band_stop_distortion_=band_stop_distortion,
    band_stop_frequency_=band_stop_frequency,
    cutoff_frequency_=cutoff_frequency)
@validate_output(butterworth_filter_order)
def calculate_butterworth_filter_order(bandwidth_distortion_: float, band_stop_distortion_: float, band_stop_frequency_: Quantity,
    cutoff_frequency_: Quantity) -> int:
    result_expr = solve(law, butterworth_filter_order, dict=True)[0][butterworth_filter_order]
    result_expr = result_expr.subs({
        bandwidth_distortion: bandwidth_distortion_,
        band_stop_distortion: band_stop_distortion_,
        band_stop_frequency: band_stop_frequency_,
        cutoff_frequency: cutoff_frequency_,
    })
    return ceil(convert_to_float(result_expr))
