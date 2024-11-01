"""
Mass flow rate
==============

Mass flow rate is the rate of change in the mass of an object. Examples include the outflow of a substance
from a certain volume, the flow in a pipe section, the combustion of fuel.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

time = symbols.time
"""
:symbols:`time`.
"""

mass_flow_rate = clone_as_function(symbols.mass_flow_rate, [time])
"""
:symbols:`mass_flow_rate` as a function of time.
"""

mass = clone_as_function(symbols.mass, [time])
"""
:symbols:`mass` as a function of time.
"""

definition = Eq(mass_flow_rate(time), Derivative(mass(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_start_=mass, mass_end_=mass, time_=time)
@validate_output(mass_flow_rate)
def calculate_mass_flow_rate(mass_start_: Quantity, mass_end_: Quantity,
    time_: Quantity) -> Quantity:
    mass_function_ = time * (mass_end_ - mass_start_) / time_
    applied_definition = definition.subs(mass(time), mass_function_)
    # calculate mass flow rate
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
