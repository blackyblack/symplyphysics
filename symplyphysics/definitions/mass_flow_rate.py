r"""
Mass flow rate
==============

Mass flow rate is the rate of change in the mass of an object. Examples include the outflow of a substance
from a certain volume, the flow in a pipe section, the combustion of fuel.

**Notation:**

#. :math:`\frac{d}{d t}` denotes a derivative w.r.t. time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input,
    validate_output)

# Definition: mu = dm / dt
# Where:
## mu is mass flow rate
## m is mass of matter (function m(t))
## t is time

mass_flow_rate = Function("mass_flow_rate", units.mass / units.time)
r"""
Mass flor rate as a function of time.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

mass_function = Function("mass_function", units.mass)
"""
:attr:`~symplyphysics.symbols.basic.mass` as a function of time.

Symbol:
    :code:`m`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(mass_flow_rate(time), Derivative(mass_function(time), time))
r"""
:code:`mu = dm/dt`

Latex:
    .. math::
        \mu = \frac{d m}{d t}
"""


@validate_input(mass_start_=units.mass, mass_end_=units.mass, time_=time)
@validate_output(mass_flow_rate)
def calculate_mass_flow_rate(mass_start_: Quantity, mass_end_: Quantity,
    time_: Quantity) -> Quantity:
    mass_function_ = time * (mass_end_ - mass_start_) / time_
    applied_definition = definition.subs(mass_function(time), mass_function_)
    # calculate mass flow rate
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
