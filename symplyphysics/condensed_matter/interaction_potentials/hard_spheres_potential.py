"""
Hard spheres potential
======================

Hard spheres are widely used as model particles in the statistical mechanical theory, defined as
impenetrable spheres that cannot overlap in space, which mimics extremely strong repulsion that
atoms and molecules experience at very close distances.

**Conditions:**

#. Spheres are identical.

**Links**:

#. `Wikipedia <https://en.wikipedia.org/wiki/Hard_spheres#Formal_definition>`__.
"""

from sympy import Eq, Piecewise, S
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

potential = symbols.potential_energy
"""
:symbols:`potential_energy` of the configuration.
"""

distance = clone_as_symbol(symbols.euclidean_distance, display_symbol="r", display_latex="r")
"""
:symbols:`euclidean_distance` between the centers of the particles.
"""

sphere_diameter = clone_as_symbol(symbols.diameter, display_symbol="sigma", display_latex="\\sigma")
"""
:symbols:`diameter` of the spheres.
"""

law = Eq(
    potential,
    Piecewise((S.Infinity, distance <= sphere_diameter), (0, distance > sphere_diameter)),
)
r"""
..
    The printers cannot yet work with `Piecewise` expressions.

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
