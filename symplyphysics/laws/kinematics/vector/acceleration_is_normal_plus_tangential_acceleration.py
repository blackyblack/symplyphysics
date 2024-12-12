"""
Acceleration is normal plus tangential acceleration
===================================================

The acceleration of a body moving arbitrarily is composed of two parts:

#. *normal, or centripetal, acceleration*, which is always present in a rotating environment
   and points to the instantaneous axis of rotation,
#. and *tangential acceleration*, which is responsible for the change in the magnitude of
   the velocity vector.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Circular_motion#Formula>`__.

#. `Mathematica LibreTexts <https://math.libretexts.org/Bookshelves/Calculus/Supplemental_Modules_(Calculus)/Vector_Calculus/2%3A_Vector-Valued_Functions_and_Motion_in_Space/2.6%3A_Tangential_and_Normal_Components_of_Acceleration>`__.
"""

from pytest import approx
from symplyphysics import (
    Quantity,
    dot_vectors,
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
    scale_vector,
)


def acceleration_law(
    normal_acceleration_: Vector,
    tangential_acceleration_: Vector,
) -> Vector:
    r"""
    Total acceleration via normal and tangential accelerations.

    Law:
        :code:`a = a_n + a_t`
    
    Latex:
        .. math::
            \vec a = {\vec a}_n + {\vec a}_\tau

    :param normal_acceleration\_: vector of normal acceleration

        Symbol: :code:`a_n`

        Latex: :math:`{\vec a}_n`

        Dimension: *acceleration*

    :param tangential_acceleration\_: vector of tangential acceleration

        Symbol: :code:`a_t`

        Latex: :math:`{\vec a}_\tau`

        Dimension: *acceleration*

    :return: vector of total acceleration

        Symbol: :code:`a`

        Latex: :math:`\vec a`

        Dimension: *acceleration*
    """

    return add_cartesian_vectors(normal_acceleration_, tangential_acceleration_)


def normal_acceleration_law(
    total_acceleration_: Vector,
    tangential_acceleration_: Vector,
) -> Vector:
    r"""
    Normal acceleration via total and tangential accelerations.

    Law:
        :code:`a_n = a - a_t`

    Latex:
        .. math::
            {\vec a}_n = \vec a - {\vec a}_\tau

    :param total_acceleration\_: vector of total acceleration

        Symbol: :code:`a`

        Latex: :math:`\vec a`

        Dimension: *acceleration*

    :param tangential_acceleration\_: vector of tangential acceleration

        Symbol: :code:`a_t`

        Latex: :math:`{\vec a}_\tau`

        Dimension: *acceleration*
    """

    opposite_tangential_vector = scale_vector(-1, tangential_acceleration_)
    return add_cartesian_vectors(total_acceleration_, opposite_tangential_vector)


def tangential_acceleration_law(
    total_acceleration_: Vector,
    normal_acceleration_: Vector,
) -> Vector:
    r"""
    Tangential acceleration via total and normal accelerations.

    Law:
        :code:`a_t = a - a_n`

    Latex:
        .. math::
            {\vec a}_\tau = \vec a - {\vec a}_n

    :param total_acceleration\_: vector of total acceleration

        Symbol: :code:`a`

        Latex: :math:`\vec a`

        Dimension: *acceleration*

    :param normal_acceleration\_: vector of normal acceleration

        Symbol: :code:`a_n`

        Latex: :math:`{\vec a}_n`

        Dimension: *acceleration*
    """

    opposite_radial_vector = scale_vector(-1, normal_acceleration_)
    return add_cartesian_vectors(total_acceleration_, opposite_radial_vector)


@validate_input(
    normal_acceleration_=units.acceleration,
    tangential_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_acceleration(
    normal_acceleration_: QuantityVector,
    tangential_acceleration_: QuantityVector,
) -> QuantityVector:
    radial_acceleration_vector = normal_acceleration_.to_base_vector()
    tangential_acceleration_vector = tangential_acceleration_.to_base_vector()
    dot_vectors_result = Quantity(
        dot_vectors(radial_acceleration_vector, tangential_acceleration_vector))
    if dot_vectors_result.scale_factor != approx(0.0, rel=1e-3):
        raise ValueError(
            "Radial and tangential acceleration vectors should be perpendicular to each other")
    acceleration_vector = acceleration_law(
        radial_acceleration_vector,
        tangential_acceleration_vector,
    )
    return QuantityVector.from_base_vector(acceleration_vector)


@validate_input(
    total_acceleration_=units.acceleration,
    tangential_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_radial_acceleration(
    total_acceleration_: QuantityVector,
    tangential_acceleration_: QuantityVector,
) -> QuantityVector:
    total_acceleration_vector = total_acceleration_.to_base_vector()
    tangential_acceleration_vector = tangential_acceleration_.to_base_vector()
    radial_acceleration_vector = normal_acceleration_law(
        total_acceleration_vector,
        tangential_acceleration_vector,
    )
    return QuantityVector.from_base_vector(radial_acceleration_vector)


@validate_input(
    total_acceleration_=units.acceleration,
    normal_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_tangential_acceleration(
    total_acceleration_: QuantityVector,
    normal_acceleration_: QuantityVector,
) -> QuantityVector:
    total_acceleration_vector = total_acceleration_.to_base_vector()
    radial_acceleration_vector = normal_acceleration_.to_base_vector()
    tangential_acceleration_vector = tangential_acceleration_law(
        total_acceleration_vector,
        radial_acceleration_vector,
    )
    return QuantityVector.from_base_vector(tangential_acceleration_vector)
