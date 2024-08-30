"""
Electric field is force over test charge
========================================

The strength of the electric field at any point in space can be found by placing a
so called test charge in that point and measuring the electrostatic force applied to
the test charge. The resulting electric field is the ratio of the force applied to
the value of the test charge.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

electric_field_strength = Symbol("electric_field_strength", units.force / units.charge)
"""
Strength of the electric field.

Symbol:
    :code:`E`
"""

electrostatic_force = clone_symbol(symbols.dynamics.force, "electrostatic_force")
"""
Projection of the electrostatic force applied to the test charge.

Symbol:
    :code:`F`
"""

test_charge = Symbol("test_charge", units.charge)
r"""
Value of the test charge.

Symbol:
    :code:`q_0`

Latex:
    :math:`q_0`
"""

law = Eq(electric_field_strength, electrostatic_force / test_charge)
r"""
:code:`E = F / q_0`

Latex:
    .. math::
        E = \frac{F}{q_0}
"""


@validate_input(electrostatic_force_=electrostatic_force, test_charge_=test_charge)
@validate_output(electric_field_strength)
def calculate_electric_field(electrostatic_force_: Quantity, test_charge_: Quantity) -> Quantity:
    result = solve(law, electric_field_strength)[0]
    result_field = result.subs({
        electrostatic_force: electrostatic_force_,
        test_charge: test_charge_,
    })
    return Quantity(result_field)
