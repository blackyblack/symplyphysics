"""
Admittance is conductance and susceptance
=========================================

Admittance is generally a complex quantity whose real part is called conductance
and whose imaginary part is called susceptance.
"""

from sympy import (Eq, I)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

admittance = symbols.admittance
"""
Admittance.
"""

conductance = symbols.conductance
"""
Conductance.
"""

susceptance = symbols.susceptance
"""
Susceptance.
"""

law = Eq(admittance, conductance + I * susceptance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    conductance_=conductance,
    susceptance_=susceptance,
)
@validate_output(admittance)
def calculate_admittance(
    conductance_: Quantity,
    susceptance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        conductance: conductance_,
        susceptance: susceptance_,
    })
    return Quantity(result)
