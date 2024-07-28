"""
Admittance is inversed impedance
================================

*Admittance, or complex conductance*, is a physical quantity measuring the
ability of a dipole to conduct electrical signal.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

dipole_admittance = Symbol("dipole_admittance", units.conductance)
"""
Admittance of the dipole.

Symbol:
    :code:`Y`
"""

dipole_impedance = Symbol("dipole_impedance", units.impedance)
"""
:doc:`Impedance <definitions.impedance_is_resistance_and_reactance>` of the dipole.

Symbol:
    :code:`Z`
"""

definition = Eq(dipole_admittance, 1 / dipole_impedance)
r"""
:code:`Y = 1 / Z`

Latex:
    .. math::
        Y = \frac{1}{Z}
"""


@validate_input(impedance_=dipole_impedance)
@validate_output(dipole_admittance)
def calculate_admittance(impedance_: Quantity) -> Quantity:
    solved = solve(definition, dipole_admittance, dict=True)[0][dipole_admittance]
    result_expr = solved.subs({dipole_impedance: impedance_})
    return Quantity(result_expr)
