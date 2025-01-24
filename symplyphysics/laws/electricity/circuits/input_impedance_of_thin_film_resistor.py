"""
Input impedance of thin film resistor
=====================================

Thin-film resistors in integrated design are used in microwave circuits. The input
impedance of a thin-film resistor depends on its resistance and capacitance, as well as
on the frequency.

..
    TODO: find link
    NOTE: replace `2 * pi * f` with `omega`
"""

from sympy import Eq, solve, I, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols

input_impedance = symbols.electrical_impedance
"""
Input :symbols:`electrical_impedance` of the resistor.
"""

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the film.
"""

frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of the current.
"""

capacitance = symbols.capacitance
"""
:symbols:`capacitance` of the film.
"""

law = Eq(input_impedance, resistance / (1 + I * 2 * pi * frequency * resistance * capacitance / 3))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    resistance_=resistance,
    frequency_=frequency,
    capacitance_=capacitance,
)
@validate_output(input_impedance)
def calculate_input_impedance(
    resistance_: Quantity,
    frequency_: Quantity,
    capacitance_: Quantity,
) -> Quantity:
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        resistance: resistance_,
        frequency: frequency_,
        capacitance: capacitance_,
    })
    return Quantity(result_expr)
