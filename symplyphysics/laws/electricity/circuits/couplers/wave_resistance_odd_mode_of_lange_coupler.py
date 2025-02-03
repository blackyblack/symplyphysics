"""
Wave impedance of odd mode of Lange coupler
===========================================

The Lange coupler is based on microstrip transmission lines. When this coupler is in
operation, both even and odd modes are distributed. Knowing the coupling coefficient
between the coupler segments, the characteristic resistance of the transmission line to
which the coupler is connected, as well as the number of coupler segments, it is
possible to calculate the wave resistance for an odd mode.

.. image:: https://habrastorage.org/r/w1560/getpro/habr/upload_files/054/d02/c8d/054d02c8d91c06425ae079d34b18ce15.jpeg
    :width: 400px
    :align: center

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    clone_as_symbol,
    SymbolNew,
)

odd_mode_wave_impedance = clone_as_symbol(symbols.wave_impedance, display_symbol="eta_o", display_latex="\\eta_\\text{o}")
"""
:symbols:`wave_impedance` of the odd mode.
"""

coupling_factor = SymbolNew("coupling_factor", dimensionless)
"""
Coupling factor between coupler segments.
"""

characteristic_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="0")
"""
:symbols:`electrical_resistance` of the transmission line.
"""

segment_count = symbols.positive_number
"""
Number of segments in Lange coupler. See :symbols:`positive_number`.
"""

_first_expression = characteristic_resistance * sqrt((1 - coupling_factor) / (1 + coupling_factor))
_second_expression = sqrt(coupling_factor**2 + (1 - coupling_factor**2) * (segment_count - 1)**2)
_third_expression = (segment_count - 1) * (1 + _second_expression) / (coupling_factor + _second_expression +
    (segment_count - 1) * (1 - coupling_factor))

law = Eq(odd_mode_wave_impedance, _first_expression * _third_expression)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(coupling_factor_=coupling_factor,
    characteristic_resistance_=characteristic_resistance,
    number_segments_=segment_count)
@validate_output(odd_mode_wave_impedance)
def calculate_wave_resistance_odd_modes(coupling_factor_: float,
    characteristic_resistance_: Quantity, number_segments_: int) -> Quantity:
    result_expr = solve(law, odd_mode_wave_impedance, dict=True)[0][odd_mode_wave_impedance]
    result_expr = result_expr.subs({
        coupling_factor: coupling_factor_,
        characteristic_resistance: characteristic_resistance_,
        segment_count: number_segments_
    })
    return Quantity(result_expr)
