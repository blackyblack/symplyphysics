"""
Total energy via relativistic mass
==================================

Fundamentally the energy of an object is synonimical to its mass.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `Wikipedia, second formula in the last group of three equations <https://en.wikipedia.org/wiki/Mass_in_special_relativity#Relativistic_energy%E2%80%93momentum_equation>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    quantities,
    validate_input,
    validate_output,
    symbols,
)

relativistic_energy = symbols.energy
"""
Total, or relativistic, :symbols:`energy` of the body.
"""

relativistic_mass = symbols.mass
"""
Relativistic :symbols:`mass` of the body.
"""

law = Eq(relativistic_energy, relativistic_mass * quantities.speed_of_light**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relativistic_mass_=relativistic_mass)
@validate_output(relativistic_energy)
def calculate_rest_energy(relativistic_mass_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_energy, dict=True)[0][relativistic_energy]
    energy_applied = result_expr.subs({relativistic_mass: relativistic_mass_})
    return Quantity(energy_applied)
