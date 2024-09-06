"""
Energy stored in capacitor via capacitance and voltage
======================================================

Capacitors store energy in the electric field when they are charged.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

energy = Symbol("energy", units.energy)
"""
Energy stored in the electric field.

Symbol:
    :code:`W`
"""

capacitance = Symbol("capacitance", units.capacitance)
"""
Capacitance of the capacitor.

Symbol:
    :code:`C`
"""

voltage = Symbol("voltage", units.voltage)
"""
Voltage across the capacitor.

Symbol:
    :code:`V`
"""

law = Eq(energy, capacitance * voltage**2 / 2)
r"""
:code:`W = 1/2 * C * V^2`

Latex:
    .. math::
        W = \frac{1}{2} C V^2
"""


@validate_input(capacitance_=capacitance, voltage_=voltage)
@validate_output(energy)
def calculate_accumulated_energy(capacitance_: Quantity, voltage_: Quantity) -> Quantity:
    result_energy_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_energy_expr.subs({capacitance: capacitance_, voltage: voltage_})
    return Quantity(result_expr)
