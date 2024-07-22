"""
Kinetic energy from mass and velocity
=====================================

The *kinetic energy* of a body is the form of energy that it possesses due to its motion.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

kinetic_energy_of_body = Symbol("kinetic_energy_of_body", units.energy)
r"""
The kinetic energy of the body.

Symbol:
    K
"""


body_velocity = Symbol("body_velocity", units.velocity)
"""
The velocity of the body.

Symbol:
    v
"""

mass = symbols.basic.mass
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the body.

Symbol:
    m
"""

law = Eq(kinetic_energy_of_body, mass * body_velocity**2 / 2)
r"""
K = m * v^2 / 2

Latex:
    :math:`K = \frac{1}{2} m v^2`
"""


@validate_input(body_mass_=mass, body_velocity_=body_velocity)
@validate_output(kinetic_energy_of_body)
def calculate_kinetic_energy(body_mass_: Quantity, body_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy_of_body, dict=True)[0][kinetic_energy_of_body]
    result_expr = result_energy_expr.subs({
        symbols.basic.mass: body_mass_,
        body_velocity: body_velocity_
    })
    return Quantity(result_expr)
