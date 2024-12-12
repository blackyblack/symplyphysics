"""
Energy stored in capacitor via capacitance and voltage
======================================================

Capacitors store energy in the electric field when they are charged.

**Links:**

#. `Wikipedia, third line of equations <https://en.wikipedia.org/wiki/Capacitor#Parallel-plate_capacitor>`__.
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
:symbols:`energy` stored in the electric field.
"""

capacitance = symbols.capacitance
"""
:symbols:`capacitance` of the capacitor.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` across the capacitor.
"""

law = Eq(energy, capacitance * voltage**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(capacitance_=capacitance, voltage_=voltage)
@validate_output(energy)
def calculate_accumulated_energy(capacitance_: Quantity, voltage_: Quantity) -> Quantity:
    result_energy_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_energy_expr.subs({capacitance: capacitance_, voltage: voltage_})
    return Quantity(result_expr)
