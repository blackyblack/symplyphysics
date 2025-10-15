"""
Relative acceleration from force and acceleration due to gravity
================================================================

Suppose a reference frame :math:`S'` is fixed to a moving body :math:`A` (e.g. Earth). For some
body :math:`B` we can write a vector equation of motion relative to :math:`S'` in the gravitational
field of body :math:`A` with the rotation of body :math:`A` taken into consideration. From this, we
can gather the meaning of the acceleration due to gravity, also known as the free fall
acceleration: it is the acceleration of body :math:`B` relative to :math:`S'` in the absence of
external forces (:math:`\\vec F = 0`) in the stationary case (the velocity of body :math:`B`
relative to :math:`S'` is zero, i.e. :math:`\\vec v = 0` and :math:`{\\vec a}_\\text{Cor} = 0`).

..
    TODO: add link to source
"""

from sympy import Eq
from symplyphysics import symbols, Quantity, validate_input, validate_output, quantities

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

relative_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_rel",
    display_latex="{\\vec a}_\\text{rel}",
)
"""
Vector of relative :symbols:`acceleration` of body :math:`B` relative to :math:`S'`
"""

force = clone_as_vector_symbol(symbols.force)
"""
Vector of the net non-gravitational :symbols:`force` exerted on body :math:`B`.
"""

mass = symbols.mass
"""
:symbols:`mass` of body :math:`B`.
"""

coriolis_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_Cor",
    display_latex="{\\vec a}_\\text{Cor}",
)
"""
Vector of the Coriolis :symbols:`acceleration` of body :math:`B`.

..
    TODO: add link to vector law
"""

acceleration_due_to_gravity = clone_as_vector_symbol(quantities.acceleration_due_to_gravity)
"""
Vector of the acceleration due to gravity.
"""

law = Eq(
    relative_acceleration,
    acceleration_due_to_gravity - coriolis_acceleration + force / mass,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    acceleration_due_to_gravity_=acceleration_due_to_gravity,
    coriolis_acceleration_=coriolis_acceleration,
    force_=force,
    mass_=mass,
)
@validate_output(relative_acceleration)
def calculate_acceleration(
    acceleration_due_to_gravity_: QuantityCoordinateVector,
    coriolis_acceleration_: QuantityCoordinateVector,
    force_: QuantityCoordinateVector,
    mass_: Quantity,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        acceleration_due_to_gravity: acceleration_due_to_gravity_,
        coriolis_acceleration: coriolis_acceleration_,
        force: force_,
        mass: mass_,
    })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 374
