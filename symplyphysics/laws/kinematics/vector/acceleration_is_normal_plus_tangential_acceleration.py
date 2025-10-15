"""
Acceleration is normal plus tangential acceleration
===================================================

The acceleration of a body moving arbitrarily is composed of two parts:

#. *normal, or centripetal, acceleration*, which is always present in a rotating environment
   and points to the instantaneous axis of rotation,
#. and *tangential acceleration*, which is responsible for the change in the magnitude of
   the velocity vector.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Velocity#Polar_coordinates>`__.

#. `Mathematica LibreTexts <https://math.libretexts.org/Bookshelves/Calculus/Supplemental_Modules_(Calculus)/Vector_Calculus/2%3A_Vector-Valued_Functions_and_Motion_in_Space/2.6%3A_Tangential_and_Normal_Components_of_Acceleration>`__.
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, symbols, Quantity

from symplyphysics.core.approx import approx_equal_numbers
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

total_acceleration = clone_as_vector_symbol(symbols.acceleration)
"""
Vector of the body's total :symbols:`acceleration`.
"""

normal_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_n",
    display_latex="{\\vec a}_\\text{n}",
)
"""
Vector of the body's normal :symbols:`acceleration`.
"""

tangential_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_t",
    display_latex="{\\vec a}_\\text{t}",
)
"""
Vector of the body's tangential :symbols:`acceleration`.
"""

law = Eq(total_acceleration, normal_acceleration + tangential_acceleration)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    normal_acceleration_=normal_acceleration,
    tangential_acceleration_=tangential_acceleration,
)
@validate_output(total_acceleration)
def calculate_acceleration(
    normal_acceleration_: QuantityCoordinateVector,
    tangential_acceleration_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    dot = Quantity(VectorDot(normal_acceleration_, tangential_acceleration_))

    if not approx_equal_numbers(dot.scale_factor, 0):
        raise ValueError("The normal and tangential accelerations must be orthogonal")

    result = solve_for_vector(law, total_acceleration).subs({
        normal_acceleration: normal_acceleration_,
        tangential_acceleration: tangential_acceleration_,
    })

    return QuantityCoordinateVector.from_expr(result)


@validate_input(
    total_acceleration_=total_acceleration,
    tangential_acceleration_=tangential_acceleration,
)
@validate_output(normal_acceleration)
def calculate_radial_acceleration(
    total_acceleration_: QuantityCoordinateVector,
    tangential_acceleration_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = solve_for_vector(law, normal_acceleration).subs({
        total_acceleration: total_acceleration_,
        tangential_acceleration: tangential_acceleration_,
    })

    return QuantityCoordinateVector.from_expr(result)


@validate_input(
    total_acceleration_=total_acceleration,
    normal_acceleration_=normal_acceleration,
)
@validate_output(tangential_acceleration)
def calculate_tangential_acceleration(
    total_acceleration_: QuantityCoordinateVector,
    normal_acceleration_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = solve_for_vector(law, tangential_acceleration).subs({
        normal_acceleration: normal_acceleration_,
        total_acceleration: total_acceleration_,
    })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 452
