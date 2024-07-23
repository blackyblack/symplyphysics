"""
Capacitance from charge and voltage
===================================

The electrical *capacitance* of an object is a physical quantity that describes the capability
of the object to store energy in the form of an electric charge. It is directly proportional to
the charge accumulated in the object and inversely proportional to the voltage across it.
"""

from sympy import Eq, solve
from symplyphysics import units, Quantity, Symbol, validate_input, validate_output

capacitance = Symbol("capacitance", units.capacitance)
"""
Capacitance of the object.

Symbol:
    :code:`C`
"""

charge = Symbol("charge", units.charge)
"""
Charge accumulated in the object.

Symbol:
    :code:`q`
"""

voltage = Symbol("voltage", units.voltage)
"""
Voltage across the object.

Symbol:
    :code:`V`
"""

definition = Eq(capacitance, charge / voltage)
r"""
:code:`C = q / V`

Latex:
    .. math::
        C = \frac{q}{V}
"""


@validate_input(charge_=charge, voltage_=voltage)
@validate_output(capacitance)
def calculate_capacitance(charge_: Quantity, voltage_: Quantity) -> Quantity:
    solved = solve(definition, capacitance, dict=True)[0][capacitance]
    result_expr = solved.subs({charge: charge_, voltage: voltage_})
    return Quantity(result_expr)
