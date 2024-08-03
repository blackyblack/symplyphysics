"""
Impedance is resistance and reactance
=====================================

*Impedance* is the combination of resistance and reactance (both inductive and capacitive) and is
a complex number, containing both real and imaginary parts. The real part of impedance is
resistance, and the imaginary part is reactance.

**Notation:**

#. :math:`i` is the imaginary unit.
"""

from sympy import (I, Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)

impedance = Symbol("impedance", units.impedance)
"""
Impedance of the system.

Symbol:
    :code:`Z`
"""

resistance = Symbol("resistance", units.impedance)
"""
Resistance of the system.

Symbol:
    :code:`R`
"""

reactance = Symbol("reactance", units.impedance)
"""
Reactance of the system.

Symbol:
    :code:`X`
"""

definition = Eq(impedance, resistance + I * reactance)
"""
:code:`Z = R + i * X`

Latex:
    .. math::
        Z = R + i X
"""


@validate_input(resistance_=resistance, reactance_=reactance)
@validate_output(impedance)
def calculate_impedance_magnitude(resistance_: Quantity, reactance_: Quantity) -> Quantity:
    solved = solve(definition, impedance, dict=True)[0][impedance]
    result_expr = solved.subs({resistance: resistance_, reactance: reactance_})
    result_magnitude = abs(result_expr)
    return Quantity(result_magnitude)
