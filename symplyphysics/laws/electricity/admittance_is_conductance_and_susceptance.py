"""
Admittance is conductance and susceptance
=========================================

Admittance is generally a complex quantity whose real part is called conductance
and whose imaginary part is called susceptance.
"""

from sympy import (Eq, I)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

admittance = Symbol("admittance", units.conductance)
"""
Admittance.

Symbol:
    :code:`Y`
"""

conductance = Symbol("conductance", units.conductance)
"""
Conductance.

Symbol:
    :code:`G`
"""

susceptance = Symbol("susceptance", units.conductance)
"""
Susceptance.

Symbol:
    :code:`B`
"""

law = Eq(admittance, conductance + I * susceptance)
r"""
:code:`Y = G + i * B`

Latex:
    .. math::
        Y = G + i B
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
