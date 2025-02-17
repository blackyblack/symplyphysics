"""
Low-pass Chebyshev filter order from distortion and frequencies
===============================================================

The approximation of the power transmission coefficient of a normalized high-pass filter
is given by approximating functions of order :math:`n`. The **Chebyshev filter** is
described by the function :math:`\\left(\\cos{n \\arccos{f^{-1}}}\\right)^{-1}` of
frequency :math:`f`. In this case, a low-pass filter is considered. The low-pass filter
lets through all frequencies from :math:`0` to the set frequency.
"""

from sympy import Eq, solve, acos, ceiling
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

filter_order = Symbol("N", dimensionless)
"""
Filter order. See :symbols:`positive_number`.
"""

bandwidth_distortion = Symbol("e", dimensionless)
"""
Bandwidth distortion, which corresponds to the number of ripples in the bandwidth.
"""

band_stop_distortion = Symbol("e_1", dimensionless)
"""
Band-stop distortion, which sets the required suppression level in the filter band-stop.
"""

cutoff_frequency = clone_as_symbol(symbols.temporal_frequency, subscript="0")
"""
Cut-off frequency, i.e. :symbols:`temporal_frequency` that corresponds to the maximum
bandwidth distortion.
"""

band_stop_frequency = clone_as_symbol(symbols.temporal_frequency, subscript="1")
"""
Band-stop frequency, i.e. :symbols:`temporal_frequency` that corresponds to the maximum
band-stop distortion.
"""

law = Eq(
    filter_order,
    acos(band_stop_distortion / bandwidth_distortion) /
    acos(band_stop_frequency / cutoff_frequency))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(bandwidth_distortion_=bandwidth_distortion,
    band_stop_distortion_=band_stop_distortion,
    band_stop_frequency_=band_stop_frequency,
    cutoff_frequency_=cutoff_frequency)
@validate_output(filter_order)
def calculate_low_pass_chebyshev_filter_order(bandwidth_distortion_: float,
    band_stop_distortion_: float, band_stop_frequency_: Quantity,
    cutoff_frequency_: Quantity) -> int:
    result_expr = solve(law, filter_order,
        dict=True)[0][filter_order]
    result_expr = result_expr.subs({
        bandwidth_distortion: bandwidth_distortion_,
        band_stop_distortion: band_stop_distortion_,
        band_stop_frequency: band_stop_frequency_,
        cutoff_frequency: cutoff_frequency_,
    })
    return int(ceiling(convert_to_float(result_expr)))
