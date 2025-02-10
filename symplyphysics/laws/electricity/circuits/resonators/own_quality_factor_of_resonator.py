"""
Quality factor of resonator
===========================

If the resonator is an oscillatory circuit, that quality factor will depend on the
resistance, inductance, and oscillation frequency.

..
    TODO: fix file name
    TODO: replace `2 * pi * f` with `omega`
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

quality_factor = symbols.quality_factor
"""
:symbols:`quality_factor` of the resonator.
"""

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the oscillating circuit.
"""

inductance = symbols.inductance
"""
:symbols:`inductance` of the oscillating circuit.
"""

frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of the current.
"""

law = Eq(quality_factor, resistance / (2 * pi * frequency * inductance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_=resistance, inductance_=inductance, frequency_=frequency)
@validate_output(quality_factor)
def calculate_quality_factor(resistance_: Quantity, inductance_: Quantity,
    frequency_: Quantity) -> float:
    result_expr = solve(law, quality_factor, dict=True)[0][quality_factor]
    result_expr = result_expr.subs({
        resistance: resistance_,
        inductance: inductance_,
        frequency: frequency_,
    })
    return convert_to_float(result_expr)
