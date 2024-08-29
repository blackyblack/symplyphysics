"""
Electrostatic potential is work to bring from infinity over charge
==================================================================

Electrostatic potential is a physical quantity defined as the amount of
work needed per unit electric charge to move it from a reference point,
usually infinity, to a specific point in an electric field.

Infinity is usually used as the reference point because this would make
the potential approach zero at an infinitely remote point.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

electrostatic_potential = Symbol("electrostatic_potential", units.voltage)
"""
Electrostatic potential of a point in an electric field.

Symbol:
    :code:`V`
"""

work = Symbol("work", units.energy)
"""
Work needed to bring the charge from the infinity.

Symbol:
    :code:`A`
"""

charge = Symbol("charge", units.charge)
"""
Value of the electric charge.

Symbol:
    :code:`q`
"""

law = Eq(electrostatic_potential, work / charge)
r"""
:code:`V = A / q`

Latex:
    .. math::
        V = \frac{A}{q}
"""


@validate_input(charge_=charge, voltage_=electrostatic_potential)
@validate_output(work)
def calculate_work(charge_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, work, dict=True)[0][work]
    result_expr = result_expr.subs({
        charge: charge_,
        electrostatic_potential: voltage_,
    })
    return Quantity(result_expr)
