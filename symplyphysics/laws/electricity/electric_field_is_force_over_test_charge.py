"""
Electric field is force over test charge
========================================

The strength of the electric field at any point in space can be found by placing a
so called test charge in that point and measuring the electrostatic force applied to
the test charge. The resulting electric field is the ratio of the force applied to
the value of the test charge.

**Links:**

#. `Wikipedia, last formula in paragraph <https://en.wikipedia.org/wiki/Electric_field#Electrostatics>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

electrostatic_force = symbols.force
"""
Projection of the electrostatic :symbols:`force` applied to the test charge.
"""

test_charge = clone_as_symbol(symbols.charge, subscript="0")
"""
Value of the test :symbols:`charge`.
"""

law = Eq(electric_field_strength, electrostatic_force / test_charge)
"""
:laws:symbol::

:laws:latex::
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
