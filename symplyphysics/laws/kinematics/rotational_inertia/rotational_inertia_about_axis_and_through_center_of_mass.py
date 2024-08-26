r"""
Rotational inertia about axis and through center of mass
========================================================

The *parallel-axis theorem* relates the rotational inertia of a body about any axis to
that of the same body about a parallel axis that extends through the body's center of mass
of mass).


**Conditions:**

#. The two axes must be parallel to each other.
#. The axis used in the calculation of :math:`I_\text{com}` must pass through the body's center of mass.
"""


from sympy import Eq
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
"""
Rotational inertia about some axis.

Symbol:
    :code:`I`
"""

rotational_inertia_through_com = Symbol("rotational_inertia_through_com",
    units.mass * units.length**2)
r"""
Rotational inertia about an axis that is parallel to the given one and passes through
the center of mass.

Symbol:
    :code:`I_com`

Latex:
    :math:`I_\text{com}`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the body.

Symbol:
    :code:`m`
"""

distance_between_axes = Symbol("distance_between_axes", units.length)
"""
Perpendicular distance between the axes.

Symbol:
    :code:`h`
"""

law = Eq(
    rotational_inertia,
    rotational_inertia_through_com + mass * distance_between_axes**2,
)
r"""
:code:`I = I_com + m * h^2`

Latex:
    .. math::
        I = I_\text{com} + m h^2
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
