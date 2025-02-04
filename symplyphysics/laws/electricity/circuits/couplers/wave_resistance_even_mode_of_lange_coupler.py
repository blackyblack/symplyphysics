"""
Wave impedance of even mode of Lange coupler
============================================

The Lange coupler is based on microstrip transmission lines. When this coupler is in
operation, both even and odd modes are distributed. Knowing the coupling coefficient
between the coupler segments, the wave impedance of the odd mode, as well as the number
of coupler segments, it is possible to calculate the wave impedance for an even mode.

.. image:: https://habrastorage.org/r/w1560/getpro/habr/upload_files/054/d02/c8d/054d02c8d91c06425ae079d34b18ce15.jpeg
    :width: 400px
    :align: center

..
    TODO: find link
    TODO: rename file to mention wave *impedance*
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    clone_as_symbol,
)

even_mode_wave_impedance = clone_as_symbol(symbols.wave_impedance, display_symbol="eta_e", display_latex="\\eta_\\text{e}")
"""
:symbols:`wave_impedance` of the even mode.
"""

odd_mode_wave_impedance = clone_as_symbol(symbols.wave_impedance, display_symbol="eta_o", display_latex="\\eta_\\text{o}")
"""
:symbols:`wave_impedance` of the odd mode.
"""

coupling_factor = SymbolNew("C", dimensionless)
"""
Coupling factor between coupler segments.
"""

segment_count = symbols.positive_number
"""
Number of segments in Lange coupler. See :symbols:`positive_number`.
"""

law = Eq(
    even_mode_wave_impedance,
    odd_mode_wave_impedance * (coupling_factor + sqrt(coupling_factor**2 +
    (1 - coupling_factor**2) * (segment_count - 1)**2)) / ((segment_count - 1) *
    (1 - coupling_factor)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(coupling_factor_=coupling_factor,
    odd_mode_wave_impedance_=odd_mode_wave_impedance,
    number_segments_=segment_count)
@validate_output(even_mode_wave_impedance)
def calculate_wave_resistance_even_modes(coupling_factor_: float,
    odd_mode_wave_impedance_: Quantity, number_segments_: int) -> Quantity:
    result_expr = solve(law, even_mode_wave_impedance, dict=True)[0][even_mode_wave_impedance]
    result_expr = result_expr.subs({
        coupling_factor: coupling_factor_,
        odd_mode_wave_impedance: odd_mode_wave_impedance_,
        segment_count: number_segments_
    })
    return Quantity(result_expr)
