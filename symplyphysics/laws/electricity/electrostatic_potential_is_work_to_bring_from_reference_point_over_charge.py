"""
Electrostatic potential is work to bring from reference point over charge
=========================================================================

Electrostatic potential is a physical quantity defined as the amount of
work needed per unit electric charge to move it from a reference point,
usually infinity, to a specific point in an electric field.

Infinity is usually used as the reference point because this would make
the potential approach zero at an infinitely remote point.

**Notes:**

#. The electric potential is defined up to a constant.

**Links:**

#. `Wikipedia, fourth formula <https://en.wikipedia.org/wiki/Electric_potential#Electrostatics>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

electrostatic_potential = symbols.electric_potential
"""
Electrostatic potential of a point in an electric field. See :symbols:`electric_potential`.
"""

work = symbols.work
"""
:symbols:`work` needed to bring the charge from the reference point.
"""

charge = symbols.charge
"""
Value of the electric :symbols:`charge`.
"""

law = Eq(electrostatic_potential, work / charge)
"""
:laws:symbol::

:laws:latex::
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
