"""
Hard spheres potential
======================

Hard spheres are widely used as model particles in the statistical mechanical theory, defined as
impenetrable spheres that cannot overlap in space, which mimics extremely strong repulsion that
atoms and molecules experience at very close distances
"""

from sympy import Eq, Piecewise, S
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

potential = Symbol("potential", units.energy)
"""
Potential energy of the configuration.

Symbol:
    :code:`U`
"""

distance = Symbol("distance", units.length)
"""
Distance between the centers of the particles.

Symbol:
    :code:`r`
"""

sphere_diameter = Symbol("sphere_diameter", units.length)
r"""
Diameter of the spheres.

Symbol:
    :code:`sigma`

Latex:
    :math:`\sigma`
"""

law = Eq(
    potential,
    Piecewise((S.Infinity, distance <= sphere_diameter), (0, distance > sphere_diameter)),
)
r"""
:code:`U = Piecewise((Infinity, r <= sigma), (0, r > sigma))`

Latex:
    .. math::
        U = \begin{cases} \infty & r \le \sigma \\ 0 & r > \sigma \end{cases}
"""


@validate_input(
    distance_=distance,
    sphere_diameter_=sphere_diameter,
)
@validate_output(potential)
def calculate_potential(
    distance_: Quantity,
    sphere_diameter_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        distance: distance_,
        sphere_diameter: sphere_diameter_,
    }).doit()
    return Quantity(result)
