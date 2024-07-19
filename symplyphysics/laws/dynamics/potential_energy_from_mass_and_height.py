"""
Potential energy from mass and height
=====================================

The potential energy is a form of energy a body possesses because of its position relative to
other objects.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, symbols)

potential_energy_of_body = Symbol("potential_energy_of_body", units.energy)
"""
The potential energy of the body.
"""

mass = symbols.basic.mass
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the body.

Symbol:
    m
"""

height = Symbol("height", units.length)

law = Eq(potential_energy_of_body, mass * units.acceleration_due_to_gravity * height)
"""
E = m * g * h

Latex:
    :math:`E = m g h`
"""


@validate_input(body_mass_=symbols.basic.mass, height_=height)
@validate_output(potential_energy_of_body)
def calculate_potential_energy(body_mass_: Quantity, height_: Quantity) -> Quantity:
    result_energy_expr = solve(law, potential_energy_of_body,
        dict=True)[0][potential_energy_of_body]
    result_expr = result_energy_expr.subs({mass: body_mass_, height: height_})
    return Quantity(result_expr)
