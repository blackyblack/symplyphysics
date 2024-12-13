r"""
Rotational inertia about axis and through center of mass
========================================================

The *parallel-axis theorem* relates the rotational inertia of a body about any axis to
that of the same body about a parallel axis that extends through the body's center of mass
of mass).

**Conditions:**

#. The two axes must be parallel to each other.
#. The axis used in the calculation of :math:`I_\text{com}` must pass through the body's center of mass.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Parallel_axis_theorem#Mass_moment_of_inertia>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

rotational_inertia = symbols.rotational_inertia
"""
:symbols:`rotational_inertia` about some axis.
"""

rotational_inertia_through_com = clone_as_symbol(
    symbols.rotational_inertia,
    display_symbol="I_com",
    display_latex="I_\\text{com}",
)
"""
:symbols:`rotational_inertia` about an axis that is parallel to the given one and passes through
the center of mass.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the body.
"""

distance_between_axes = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the axes.
"""

law = Eq(
    rotational_inertia,
    rotational_inertia_through_com + mass * distance_between_axes**2,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    rotational_inertia_through_com_=rotational_inertia_through_com,
    mass_=mass,
    distance_between_axes_=distance_between_axes,
)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(
    rotational_inertia_through_com_: Quantity,
    mass_: Quantity,
    distance_between_axes_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia_through_com: rotational_inertia_through_com_,
        mass: mass_,
        distance_between_axes: distance_between_axes_,
    })
    return Quantity(result)
