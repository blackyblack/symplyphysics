"""
Falling body displacement
=========================

Suppose a reference frame :math:`S'` is fixed to a moving object :math:`A` (e.g., Earth) and some
body :math:`B` moving freely (i.e. the sum of external non-gravitational forces acting on it is
zero). In the case of an inertial frame of reference, the displacement of body :math:`B` from the
starting position would follow the usual rule :math:`\\vec s = {\\vec v}_0 t + \\frac{1}{2} {\\vec g} t^2`.
But in the case of non-inertial frames, we have to take the Coriolis force and the centrifugal force
into account as well, which results into the series shown below.

**Notes:**

#. Note that the series is truncated at the fifth term. More terms can be obtained by plugging the
   result into the equation of motion :math:`\\vec a = \\vec g + \\left[ \\vec v, \\vec \\omega \\right]`
   and integrating it over time.

**Conditions:**

#. The sum :math:`\\vec F` of all other, non-gravitational forces acting on body :math:`B` is
   :math:`0`.
#. :math:`\\vec g` is independent of coordinates.

..
    TODO: add link to source
    TODO: add `O(t^5)` to the law
"""

from sympy import Eq, evaluate
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

displacement = clone_as_vector_symbol(symbols.distance)
"""
Vector of the displacement of body :math:`B` in frame :math:`S'`. See :symbols:`distance`.
"""

time = symbols.time
"""
:symbols:`time`.
"""

initial_velocity = clone_as_vector_symbol(symbols.speed, subscript="0")
"""
Vector of the initial velocity of body :math:`B`. See :symbols:`speed`.
"""

angular_velocity = clone_as_vector_symbol(symbols.angular_speed)
"""
Pseudovector of the angular velocity of body :math:`B`. See :symbols:`angular_speed`.
"""

acceleration_due_to_gravity = clone_as_vector_symbol(quantities.acceleration_due_to_gravity)
"""
Vector of the acceleration due to gravity of body :math:`B`.
"""

with evaluate(False):
    _first_power = initial_velocity * time

    _second_power = time**2 * (acceleration_due_to_gravity / 2 +
        VectorCross(initial_velocity, angular_velocity))

    _third_power = (time**3 / 3) * (VectorCross(acceleration_due_to_gravity, angular_velocity) +
        2 * VectorCross(VectorCross(initial_velocity, angular_velocity), angular_velocity))

    _fourth_power = (time**4 / 6) * VectorCross(
        VectorCross(acceleration_due_to_gravity, angular_velocity), angular_velocity)

law = Eq(
    displacement,
    _first_power + _second_power + _third_power + _fourth_power,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    time_=time,
    initial_velocity_=initial_velocity,
    angular_velocity_=angular_velocity,
    acceleration_due_to_gravity_=acceleration_due_to_gravity,
)
@validate_output(displacement)
def calculate_displacement(
    time_: Quantity,
    initial_velocity_: QuantityCoordinateVector,
    angular_velocity_: QuantityCoordinateVector,
    acceleration_due_to_gravity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        time: time_,
        initial_velocity: initial_velocity_,
        angular_velocity: angular_velocity_,
        acceleration_due_to_gravity: acceleration_due_to_gravity_,
    })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 372
