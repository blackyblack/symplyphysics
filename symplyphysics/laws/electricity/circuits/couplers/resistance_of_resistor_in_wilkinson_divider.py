"""
Resistor resistance in Wilkinson divider
========================================

The Wilkinson divider is a device designed to divide the power of a microwave signal
into two output ports. The microstrip version has one surface-mounted resistor. The
resistance of this resistor depends on the resistance of the transmission line and the
power ratio at the output ports.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    clone_as_symbol,
)

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the resistor.
"""

transmission_line_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="0")
"""
:symbols:`electrical_resistance` of the transmission line.
"""

power_ratio = SymbolNew("k", dimensionless)
"""
Ratio of the power at the outputs of the divider.
"""

law = Eq(resistance, transmission_line_resistance * (1 + power_ratio**2) / power_ratio)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(transmission_line_resistance_=transmission_line_resistance,
    ratio_of_power_=power_ratio)
@validate_output(resistance)
def calculate_resistance(transmission_line_resistance_: Quantity,
    ratio_of_power_: float) -> Quantity:
    result_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_expr.subs({
        transmission_line_resistance: transmission_line_resistance_,
        power_ratio: ratio_of_power_,
    })
    return Quantity(result_expr)
