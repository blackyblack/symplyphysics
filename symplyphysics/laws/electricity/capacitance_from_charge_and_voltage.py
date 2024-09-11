"""
Capacitance from charge and voltage
===================================

The electrical *capacitance* of an object is a physical quantity that describes the capability
of the object to store energy in the form of an electric charge. It is directly proportional to
the charge accumulated in the object and inversely proportional to the voltage across it.
"""

from sympy import Eq, solve
from symplyphysics import units, Quantity, SymbolNew, validate_input, validate_output

capacitance = SymbolNew("C", units.capacitance)
"""
Capacitance of the object.
"""

charge = SymbolNew("q", units.charge)
"""
Charge accumulated in the object.
"""

voltage = SymbolNew("V", units.voltage)
"""
Voltage across the object.
"""

definition = Eq(capacitance, charge / voltage)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_=charge, voltage_=voltage)
@validate_output(capacitance)
def calculate_capacitance(charge_: Quantity, voltage_: Quantity) -> Quantity:
    solved = solve(definition, capacitance, dict=True)[0][capacitance]
    result_expr = solved.subs({charge: charge_, voltage: voltage_})
    return Quantity(result_expr)
