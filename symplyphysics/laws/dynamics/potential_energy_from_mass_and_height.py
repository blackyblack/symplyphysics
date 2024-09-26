"""
Potential energy from mass and height
=====================================

The potential energy is a form of energy a body possesses because of its position relative to
other objects.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

potential_energy = symbols.potential_energy
"""
The :symbols:`potential_energy` of the body.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the body.
"""

height = symbols.height
"""
The elevation from ground level. See :symbols:`height`
"""

law = Eq(potential_energy, mass * quantities.acceleration_due_to_gravity * height)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(body_mass_=mass, height_=height)
@validate_output(potential_energy)
def calculate_potential_energy(body_mass_: Quantity, height_: Quantity) -> Quantity:
    result_energy_expr = solve(law, potential_energy, dict=True)[0][potential_energy]
    result_expr = result_energy_expr.subs({mass: body_mass_, height: height_})
    return Quantity(result_expr)
