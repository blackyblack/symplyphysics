"""
High pass Chebyshev filter from distortion and frequency
========================================================

The approximation of the power transmission coefficient of a normalized high-pass filter
is given by approximating functions of order :math:`n`. The **Chebyshev filter** is
described by the function :math:`\\left(\\cos{n \\acos{f^{-1}}}\\right)^{-1}` of
frequency :math:`f`. In this case, a high-pass filter is considered. The high-pass
filter does removes all frequencies from :math:`0` to the set frequency.

..
    TODO: find link
"""

from sympy import Eq, solve, acosh, ceiling
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    clone_as_symbol,
    symbols,
)

filter_order = SymbolNew("N", dimensionless)
"""
Filter order. See :symbols:`positive_number`.
"""

bandwidth_distortion = SymbolNew("e", dimensionless)
"""
Bandwidth distortion, which corresponds to the number of ripples in the bandwidth.
"""

band_stop_distortion = SymbolNew("e_1", dimensionless)
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
    acosh(band_stop_distortion / bandwidth_distortion) /
    acosh(cutoff_frequency / band_stop_frequency))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(bandwidth_distortion_=bandwidth_distortion,
    band_stop_distortion_=band_stop_distortion,
    band_stop_frequency_=band_stop_frequency,
    cutoff_frequency_=cutoff_frequency)
@validate_output(filter_order)
def calculate_chebyshev_filter_order(bandwidth_distortion_: float, band_stop_distortion_: float,
    band_stop_frequency_: Quantity, cutoff_frequency_: Quantity) -> int:
    result_expr = solve(law, filter_order, dict=True)[0][filter_order]
    result_expr = result_expr.subs({
        bandwidth_distortion: bandwidth_distortion_,
        band_stop_distortion: band_stop_distortion_,
        band_stop_frequency: band_stop_frequency_,
        cutoff_frequency: cutoff_frequency_,
    })
    return int(ceiling(convert_to_float(result_expr)))
