"""
Admittance is inversed impedance
================================

*Admittance, or complex conductance*, is a physical quantity measuring the
ability of a circuit or device to conduct electrical current.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

admittance = Symbol("admittance", units.conductance)
"""
Admittance of the object.

Symbol:
    :code:`Y`
"""

impedance = Symbol("impedance", units.impedance)
"""
:doc:`Impedance <definitions.impedance_is_resistance_and_reactance>` of the object.

Symbol:
    :code:`Z`
"""

definition = Eq(admittance, 1 / impedance)
r"""
:code:`Y = 1 / Z`

Latex:
    .. math::
        Y = \frac{1}{Z}
"""


@validate_input(impedance_=impedance)
@validate_output(admittance)
def calculate_admittance(impedance_: Quantity) -> Quantity:
    solved = solve(definition, admittance, dict=True)[0][admittance]
    result_expr = solved.subs({impedance: impedance_})
    return Quantity(result_expr)
