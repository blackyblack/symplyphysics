"""
Cross section of interaction in model of quasirigid spheres
===========================================================

The effective cross section is a physical quantity characterizing the probability of
transition of a system of two interacting particles to a certain final state, a
quantitative characteristic of the acts of collision of particles of a stream hitting a
target with target particles. The effective cross-section has the dimension of the area.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Cross_section_(physics)#Collision_among_gas_particles>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

cross_section = symbols.cross_section
"""
:symbols:`cross_section` of interaction of particles.
"""

distance_of_convergence = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` of greatest convergence of two colliding particles.
"""

law = Eq(cross_section, pi * distance_of_convergence**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(distance_of_convergence_=distance_of_convergence)
@validate_output(cross_section)
def calculate_cross_sectional_area_of_interaction(
        distance_of_convergence_: Quantity) -> Quantity:
    result_expr = solve(law, cross_section,
        dict=True)[0][cross_section]
    result_expr = result_expr.subs({
        distance_of_convergence: distance_of_convergence_,
    })
    return Quantity(result_expr)
