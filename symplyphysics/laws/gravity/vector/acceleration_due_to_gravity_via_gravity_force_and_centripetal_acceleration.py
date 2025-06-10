"""
Acceleration due to gravity via gravity force and centripetal acceleration
==========================================================================

Suppose a reference frame :math:`S'` is fixed to a rotating body :math:`A` (e.g. Earth), so that
frame :math:`S'` rotates w.r.t. another static reference frame :math:`S`. The acceleration due to
gravity (in moving frame :math:`S'`) is the acceleration another body :math:`B` has in the gravity
field of body :math:`A`, with rotational effects such as the centripetal acceleration taken into
account. It is the same for all bodies at a fixed point, but can be different at different points
in space.

..
    TODO: add link to source
"""

from sympy import Eq
from symplyphysics import symbols, Quantity, validate_input, validate_output, quantities

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

mass = symbols.mass
"""
:symbols:`mass` of body :math:`B`.
"""

acceleration_due_to_gravity = clone_as_vector_symbol(quantities.acceleration_due_to_gravity)
"""
Vector of acceleration due to gravity of body :math:`B`.
"""

gravity_force = clone_as_vector_symbol(symbols.force)
"""
Vector of the :symbols:`force` of gravity pull exerted on body :math:`B`.
"""

centripetal_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_cp",
    display_latex="{\\vec a}_\\text{cp}",
)
"""
Vector of centripetal :symbols:`acceleration` of body :math:`B`.
"""

law = Eq(acceleration_due_to_gravity, gravity_force / mass - centripetal_acceleration)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    gravity_force_=gravity_force,
    centripetal_acceleration_=centripetal_acceleration,
    mass_=mass,
)
@validate_output(acceleration_due_to_gravity)
def calculate_acceleraton_due_to_gravity(
    gravity_force_: QuantityCoordinateVector,
    centripetal_acceleration_: QuantityCoordinateVector,
    mass_: Quantity,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        mass: mass_,
        gravity_force: gravity_force_,
        centripetal_acceleration: centripetal_acceleration_,
    })

    return QuantityCoordinateVector.from_expr(result)
