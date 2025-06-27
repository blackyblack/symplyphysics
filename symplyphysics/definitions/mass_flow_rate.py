"""
Mass flow rate
==============

The *mass flow rate* of a system is the time derivative of its mass. Typical
examples include fluid passing through a pipe section or fuel being burned.

**Conditions:**

#. Mass is a continuous, differentiable function of time.

**Links:**

#. `Wikipedia – Formulation <https://en.wikipedia.org/wiki/Mass_flow_rate#Formulation>`__
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
    # mass changes linearly over the time interval `time_`
    mass_function_ = time * (mass_end_ - mass_start_) / time_
    applied_definition = definition.subs(mass(time), mass_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
