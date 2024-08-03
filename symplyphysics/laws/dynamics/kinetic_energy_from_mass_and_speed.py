"""
Kinetic energy from mass and speed
==================================

The *kinetic energy* of a body is the form of energy that it possesses due to its motion.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

kinetic_energy = Symbol("kinetic_energy", units.energy)
r"""
The kinetic energy of the body.

Symbol:
    :code:`K`
"""

speed = Symbol("speed", units.velocity)
"""
The speed of the body.

Symbol:
    :code:`v`
"""

mass = symbols.basic.mass
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the body.

Symbol:
    :code:`m`
"""

law = Eq(kinetic_energy, mass * speed**2 / 2)
r"""
:code:`K = m * v^2 / 2`

Latex:
    :math:`K = \frac{1}{2} m v^2`
"""


@validate_input(body_mass_=mass, body_velocity_=speed)
@validate_output(kinetic_energy)
def calculate_kinetic_energy(body_mass_: Quantity, body_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
    result_expr = result_energy_expr.subs({mass: body_mass_, speed: body_velocity_})
    return Quantity(result_expr)
