"""
Admittance is inverse impedance
===============================

*Admittance* (sometimes called complex conductance) describes how easily alternating current flows through a network. It is defined as the reciprocal of impedance.

Also see :doc:`Impedance law <definitions.impedance_is_resistance_and_reactance>`

**Conditions:**

#. Applicable under sinusoidal steady-state (phasor-domain) analysis for linear, time-invariant systems.

**Links:**

#. `Wikipedia – Admittance <https://en.wikipedia.org/wiki/Admittance>`__
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

admittance = symbols.admittance
"""
:symbols:`admittance` of the object.
"""

impedance = symbols.electrical_impedance
"""
:symbols:`electrical_impedance` of the object.
"""

definition = Eq(admittance, 1 / impedance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(impedance_=impedance)
@validate_output(admittance)
def calculate_admittance(impedance_: Quantity) -> Quantity:
    solved = solve(definition, admittance, dict=True)[0][admittance]
    result_expr = solved.subs({impedance: impedance_})
    return Quantity(result_expr)
