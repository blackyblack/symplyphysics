"""
Rotational inertia of a particle
================================

Rotational inertia of a rotating particle is defined as the product of its mass
and the squared radius of rotation.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
"""
Rotational inertia of the particle.

Symbol:
    :code:`I`
"""

distance_to_axis = Symbol("distance_to_axis", units.length)
"""
Distance to axis, or radius of rotation.

Symbol:
    :code:`r`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the particle.
"""

law = Eq(rotational_inertia, mass * distance_to_axis**2)
r"""
:code:`I = m * r^2`

Latex:
    .. math::
        I = m r^2
"""


@validate_input(mass_=mass, radius_=distance_to_axis)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        distance_to_axis: radius_,
    })
    return Quantity(result)
