"""
Input impedance of thin film resistor
=====================================

Thin-film resistors in integrated design are used in microwave circuits. The input
impedance of a thin-film resistor depends on its resistance and capacitance, as well as
on the angular frequency.

..
    TODO: find link
"""

from sympy import Eq, solve, I
from symplyphysics import Quantity, validate_input, validate_output, symbols

input_impedance = symbols.electrical_impedance
"""
Input :symbols:`electrical_impedance` of the resistor.
"""

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the film.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the current.
"""

capacitance = symbols.capacitance
"""
:symbols:`capacitance` of the film.
"""

law = Eq(input_impedance, resistance / (1 + I * angular_frequency * resistance * capacitance / 3))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    resistance_=resistance,
    frequency_=angular_frequency,
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
        angular_frequency: frequency_,
        capacitance: capacitance_,
    })
    return Quantity(result_expr)
