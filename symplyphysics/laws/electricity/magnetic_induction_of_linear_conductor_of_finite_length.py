"""
Magnetic flux density of linear conductor of finite length
==========================================================

Let there be a rectilinear conductor of finite length. Then its magnetic flux density
will depend on the magnitude of the current and the material. It also depends on the
perpendicular distance to the conductor and on the angles between the lines drawn from
the ends of the conductor to the point and the conductor.

**Conditions:**

#. Conductor should be rectilinear.
#. Length of the conductor is finite.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, cos
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

magnetic_flux_density = symbols.magnetic_flux_density
"""
:symbols:`magnetic_flux_density` through the conductor.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the medium.
"""

current = symbols.current
"""
:symbols:`current` running through the conductor.
"""

first_angle = clone_as_symbol(symbols.angle, subscript="1")
"""
:symbols:`angle` between origin and the first end of the conductor.
"""

second_angle = clone_as_symbol(symbols.angle, subscript="2")
"""
:symbols:`angle` between origin and the second end of the conductor.
"""

distance = symbols.euclidean_distance
"""
Perpendicular :symbols:`euclidean_distance` to the conductor.
"""

law = Eq(
    magnetic_flux_density,
    absolute_permeability * current * (cos(first_angle) + cos(second_angle)) /
    (4 * pi * distance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permeability_=absolute_permeability,
    current_=current,
    first_angle_=first_angle,
    second_angle_=second_angle,
    distance_=distance)
@validate_output(magnetic_flux_density)
def calculate_induction(relative_permeability_: float, current_: Quantity,
    first_angle_: float | Quantity, second_angle_: float | Quantity,
    distance_: Quantity) -> Quantity:
    result_expr = solve(law, magnetic_flux_density, dict=True)[0][magnetic_flux_density]
    result_expr = result_expr.subs({
        absolute_permeability: relative_permeability_,
        current: current_,
        first_angle: first_angle_,
        second_angle: second_angle_,
        distance: distance_
    })
    return Quantity(result_expr)
