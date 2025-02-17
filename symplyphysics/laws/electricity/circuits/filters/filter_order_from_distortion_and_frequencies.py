"""
Filter order from distortion and frequency
==========================================

The approximation of the power transmission coefficient of a normalized low-pass filter
is given by approximating functions of the order of :math:`n`.

..
    TODO: find link
"""

from sympy import Expr, Eq, solve, ceiling
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

filter_function = Symbol("F", dimensionless)
"""
Approximating function of order :math:`n` (:attr:`~filter_order`) that depends on
:symbols:`temporal_frequency`.
"""

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

band_stop_frequency = clone_as_symbol(symbols.temporal_frequency, subscript="1")
"""
Band-stop frequency, i.e. :symbols:`temporal_frequency` that corresponds to the maximum
band-stop distortion.
"""

law = Eq(filter_function, band_stop_distortion / bandwidth_distortion)
"""
:laws:symbol::

:laws:latex::
"""

@validate_input(band_stop_frequency_=band_stop_frequency,
    bandwidth_distortion_=bandwidth_distortion,
    band_stop_distortion_=band_stop_distortion)
@validate_output(filter_order)
def calculate_order(filter_function_: Expr, band_stop_frequency_: Quantity,
    bandwidth_distortion_: float, band_stop_distortion_: float) -> int:
    result_expr = law.subs({
        filter_function: filter_function_,
        band_stop_frequency: band_stop_frequency_,
        bandwidth_distortion: bandwidth_distortion_,
        band_stop_distortion: band_stop_distortion_,
    })
    result = solve(result_expr, filter_order, dict=True)[0][filter_order]
    return int(ceiling(convert_to_float(result)))
