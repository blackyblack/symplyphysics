"""
Band pass Chebyshev filter order from distortion and frequency
==============================================================

The approximation of the power transmission coefficient of a normalized high-pass filter
is given by approximating functions of order :math:`n`. The **Chebyshev filter** is
described by the function :math:`\\left(\\cos{n \\acos{f^{-1}}}\\right)^{-1}` of
frequency :math:`f`. In this case, a band-pass filter is considered. The bandpass filter
only passes frequencies from a certain frequency range.

**Notes:**

#. Also refer to `this link (in Russian) <http://www.dsplib.ru/content/filters/ch8/ch8.html#r4>`__.

..
    TODO: find link
    TODO: fix file name
"""

from sympy import Eq, solve, acosh, ceiling
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

filter_order = SymbolNew("N", dimensionless)
"""
Chebyshev filter order. See :symbols:`positive_number`.
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

bandwidth = clone_as_symbol(symbols.temporal_frequency, display_symbol="Delta(f)", display_latex="\\Delta f")
"""
Bandwidth. See :symbols:`temporal_frequency`.
"""

law = Eq(
    filter_order,
    acosh(bandwidth_distortion / band_stop_distortion) / acosh(
    (band_stop_frequency**2 - cutoff_frequency**2) / (bandwidth * band_stop_frequency)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(bandwidth_distortion_=bandwidth_distortion,
    band_stop_distortion_=band_stop_distortion,
    band_stop_frequency_=band_stop_frequency,
    cutoff_frequency_=cutoff_frequency,
    bandwidth_=bandwidth)
@validate_output(filter_order)
def calculate_band_pass_chebyshev_filter_order(bandwidth_distortion_: float,
    band_stop_distortion_: float, band_stop_frequency_: Quantity, cutoff_frequency_: Quantity,
    bandwidth_: Quantity) -> int:
    result_expr = solve(law, filter_order,
        dict=True)[0][filter_order]
    result_expr = result_expr.subs({
        bandwidth_distortion: bandwidth_distortion_,
        band_stop_distortion: band_stop_distortion_,
        band_stop_frequency: band_stop_frequency_,
        cutoff_frequency: cutoff_frequency_,
        bandwidth: bandwidth_,
    })
    return int(ceiling(convert_to_float(result_expr)))
