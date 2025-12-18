"""
Rotational inertia of a particle
================================

Rotational inertia of a rotating particle is defined as the product of its mass
and the squared radius of rotation.

**Links:**

#. `Wikipedia, third equation <https://en.wikipedia.org/wiki/Moment_of_inertia#Definition>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

rotational_inertia = symbols.rotational_inertia
"""
:symbols:`rotational_inertia` of the particle.
"""

distance_to_axis = symbols.distance_to_axis
"""
:symbols:`distance_to_axis`, or radius of rotation.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the particle.
"""

law = Eq(rotational_inertia, mass * distance_to_axis**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_=mass, radius_=distance_to_axis)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        distance_to_axis: radius_,
    })
    return Quantity(result)
