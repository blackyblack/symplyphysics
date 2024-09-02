"""
Potential energy from mass and height
=====================================

The potential energy is a form of energy a body possesses because of its position relative to
other objects.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

potential_energy = Symbol("potential_energy", units.energy)
"""
The potential energy of the body.

Symbol:
    :code:`U`
"""

mass = symbols.basic.mass
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the body.
"""

height = Symbol("height", units.length)
"""
The elevation from ground level.

Symbol:
    :code:`h`
"""

law = Eq(potential_energy, mass * units.acceleration_due_to_gravity * height)
"""
:code:`U = m * g * h`

Latex:
    .. math::
        U = m g h
"""


@validate_input(body_mass_=mass, height_=height)
@validate_output(potential_energy)
def calculate_potential_energy(body_mass_: Quantity, height_: Quantity) -> Quantity:
    result_energy_expr = solve(law, potential_energy, dict=True)[0][potential_energy]
    result_expr = result_energy_expr.subs({mass: body_mass_, height: height_})
    return Quantity(result_expr)
