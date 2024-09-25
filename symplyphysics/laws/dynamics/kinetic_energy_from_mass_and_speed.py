"""
Kinetic energy from mass and speed
==================================

The *kinetic energy* of a body is the form of energy that it possesses due to its motion.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

kinetic_energy = symbols.kinetic_energy
r"""
The :symbols:`kinetic_energy` of the body.
"""

speed = symbols.speed
"""
The :symbols:`speed` of the body.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the body.
"""

law = Eq(kinetic_energy, mass * speed**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(body_mass_=mass, body_velocity_=speed)
@validate_output(kinetic_energy)
def calculate_kinetic_energy(body_mass_: Quantity, body_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
    result_expr = result_energy_expr.subs({mass: body_mass_, speed: body_velocity_})
    return Quantity(result_expr)
