"""
Mass flow rate
==============

Mass flow rate is the rate of change in the mass of an object. Examples include the outflow of a substance
from a certain volume, the flow in a pipe section, the combustion of fuel.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, FunctionNew, SymbolNew, validate_input, validate_output)

mass_flow_rate = FunctionNew("mu(t)", units.mass / units.time, display_latex="\\mu")
"""
Mass flow rate as a function of time.
"""

mass = FunctionNew("m(t)", units.mass, display_latex="m")
"""
:attr:`~symplyphysics.symbols.mass` as a function of time.
"""

time = SymbolNew("t", units.time)
"""
Time.
"""

definition = Eq(mass_flow_rate(time), Derivative(mass(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_start_=units.mass, mass_end_=units.mass, time_=time)
@validate_output(mass_flow_rate)
def calculate_mass_flow_rate(mass_start_: Quantity, mass_end_: Quantity,
    time_: Quantity) -> Quantity:
    mass_function_ = time * (mass_end_ - mass_start_) / time_
    applied_definition = definition.subs(mass(time), mass_function_)
    # calculate mass flow rate
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
