"""
Standing wave ratio from voltage
================================

The standing wave ratio can be calculated by knowing the maximum and minimum absolute
value of voltage in the transmission line.

**Links:**

#. `Engineering LibreTexts <https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Electro-Optics/Book%3A_Electromagnetics_I_(Ellingson)/03%3A_Transmission_Lines/3.14%3A_Standing_Wave_Ratio>`__.

..
    TODO: fix file name
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

standing_wave_ratio = symbols.standing_wave_ratio
"""
:symbols:`standing_wave_ratio`.
"""

maximum_voltage_module = clone_as_symbol(symbols.voltage, display_symbol="min(abs(V))", display_latex="\\min{|V|}")
"""
Minimum absolute value of :symbols:`voltage` in the transmission line.
"""

minimum_voltage_module = clone_as_symbol(symbols.voltage, display_symbol="max(abs(V))", display_latex="\\max{|V|}")
"""
Maximum absolute value of :symbols:`voltage` in the transmission line.
"""

law = Eq(standing_wave_ratio, maximum_voltage_module / minimum_voltage_module)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(maximum_voltage_=maximum_voltage_module, minimum_voltage_=minimum_voltage_module)
@validate_output(standing_wave_ratio)
def calculate_coefficient_standing_wave(maximum_voltage_: Quantity,
    minimum_voltage_: Quantity) -> float:
    if maximum_voltage_.scale_factor < minimum_voltage_.scale_factor:
        raise ValueError("The maximum voltage must be not less than the minimum voltage")
    result_expr = solve(law, standing_wave_ratio, dict=True)[0][standing_wave_ratio]
    result_expr = result_expr.subs({
        maximum_voltage_module: maximum_voltage_,
        minimum_voltage_module: minimum_voltage_,
    })
    return convert_to_float(result_expr)
