"""
Energy stored in inductor via inductance and current
====================================================

Inductors store energy in the magnetic field when the current is flowing
through it.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

energy = symbols.work
"""
Energy stored in the magnetic field.
"""

inductance = symbols.inductance
"""
Inductance of the inductor.
"""

current = symbols.current
"""
Current flowing through the inductor.
"""

law = Eq(energy, inductance * current**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(inductance_=inductance, current_=current)
@validate_output(energy)
def calculate_accumulated_energy(inductance_: Quantity, current_: Quantity) -> Quantity:
    result_energy_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_energy_expr.subs({inductance: inductance_, current: current_})
    return Quantity(result_expr)
