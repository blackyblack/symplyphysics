"""
Time constant of resistor-capacitor circuit
===========================================

The time constant of an RC circuit is the time it takes for the current to become
:math:`e` times smaller its original value. Also see `exponential decay
<https://en.wikipedia.org/wiki/Exponential_decay>`_.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

time_constant = symbols.electric_time_constant
"""
:symbols:`electric_time_constant` of the circuit.
"""

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the resistor.
"""

capacitance = symbols.capacitance
"""
:symbols:`capacitance` of the capacitor.
"""

law = Eq(time_constant, resistance * capacitance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_=resistance, capacitance_=capacitance)
@validate_output(time_constant)
def calculate_time_constant(resistance_: Quantity, capacitance_: Quantity) -> Quantity:
    result = law.rhs.subs({
        resistance: resistance_,
        capacitance: capacitance_,
    })
    return Quantity(result)
