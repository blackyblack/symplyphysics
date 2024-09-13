"""
Coil impedance via inductance and frequency
===========================================

The impedance of a coil depends on its inductance and the angular frequency of the alternating
current. It is a purely imaginary quantity since the resistance of a coil is zero.
"""

from sympy import I, Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

impedance = symbols.electrical_impedance
"""
Impedance of the coil.
"""

angular_frequency = symbols.angular_frequency
"""
Angular frequency of the current.
"""

inductance = symbols.inductance
"""
Coil inductance.
"""

law = Eq(impedance, I * angular_frequency * inductance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(inductance_=inductance, circular_frequency_=angular_frequency)
@validate_output(impedance)
def calculate_impedance(inductance_: Quantity, circular_frequency_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, impedance, dict=True)[0][impedance]
    result_expr = result_impedance_expr.subs({
        inductance: inductance_,
        angular_frequency: circular_frequency_
    })
    return Quantity(result_expr)
