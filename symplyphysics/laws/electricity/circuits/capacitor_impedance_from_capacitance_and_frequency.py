"""
Capacitor impedance from capacitance and frequency
==================================================

The impedance of the capacitor is a purely imaginary reactive impedance which
is inversely proportional to its capacitance. Note that the resistance of the
capacitor is infinite, i.e. it is considered to be an opened connection.
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
:symbols:`electrical_impedance` of the capacitor.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the alternating current.
"""

capacitance = symbols.capacitance
"""
:symbols:`capacitance` of the capacitor.
"""

law = Eq(impedance, -I / (angular_frequency * capacitance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(capacitance_=capacitance, circular_frequency_=angular_frequency)
@validate_output(impedance)
def calculate_impedance(capacitance_: Quantity, circular_frequency_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, impedance, dict=True)[0][impedance]
    result_expr = result_impedance_expr.subs({
        capacitance: capacitance_,
        angular_frequency: circular_frequency_
    })
    return Quantity(result_expr)
