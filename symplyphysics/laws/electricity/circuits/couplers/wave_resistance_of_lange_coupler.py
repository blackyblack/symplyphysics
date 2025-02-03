"""
Wave resistance of Lange coupler
================================

The Lange coupler is based on microstrip transmission lines. When this coupler is in
operation, both even and odd modes are distributed. Knowing the wave resistance for even
and odd modes, it is possible to calculate the equivalent wave resistance of the
coupler.

.. image:: https://habrastorage.org/r/w1560/getpro/habr/upload_files/054/d02/c8d/054d02c8d91c06425ae079d34b18ce15.jpeg
    :width: 400px
    :align: center

..
    TODO: find link
    TODO: rename file
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

wave_impedance = symbols.wave_impedance
"""
:symbols:`wave_impedance` of the Lange coupler.
"""

odd_mode_wave_impedance = clone_as_symbol(symbols.wave_impedance, display_symbol="eta_o", display_latex="\\eta_\\text{o}")
"""
:symbols:`wave_impedance` of the odd mode.
"""

even_mode_wave_impedance = clone_as_symbol(symbols.wave_impedance, display_symbol="eta_e", display_latex="\\eta_\\text{e}")
"""
:symbols:`wave_impedance` of the even mode.
"""

segment_count = symbols.positive_number
"""
Number of the segments of the Lange coupler. See :symbols:`positive_number`.
"""

law = Eq(
    wave_impedance,
    sqrt((odd_mode_wave_impedance * even_mode_wave_impedance *
    (odd_mode_wave_impedance + even_mode_wave_impedance)**2) /
    ((odd_mode_wave_impedance + even_mode_wave_impedance * (segment_count - 1)) *
    (even_mode_wave_impedance + odd_mode_wave_impedance * (segment_count - 1)))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wave_resistance_odd_modes_=odd_mode_wave_impedance,
    wave_resistance_even_modes_=even_mode_wave_impedance,
    number_segments_=segment_count)
@validate_output(wave_impedance)
def calculate_wave_resistance(wave_resistance_odd_modes_: Quantity,
    wave_resistance_even_modes_: Quantity, number_segments_: int) -> Quantity:
    result_expr = solve(law, wave_impedance, dict=True)[0][wave_impedance]
    result_expr = result_expr.subs({
        odd_mode_wave_impedance: wave_resistance_odd_modes_,
        even_mode_wave_impedance: wave_resistance_even_modes_,
        segment_count: number_segments_
    })
    return Quantity(result_expr)
