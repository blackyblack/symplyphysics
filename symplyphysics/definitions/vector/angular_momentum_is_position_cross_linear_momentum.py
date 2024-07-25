r"""
Angular momentum is position cross linear momentum
==================================================

The pseudovector of *angular momentum* of a particle is the cross product of its
position vector and linear momentum. Unlike linear momentum, angular momentum depends
on where the origin of the coordinate system is chosen since it depends on the position
vector of the particle defined relative to that origin.

**Notation:**

#. :math:`\vec a \times \vec b` (:code:`cross(a, b)`) denotes a cross product between
   vectors :math:`\vec a` and :math:`\vec b`.
"""

from symplyphysics import (
    units,
    Vector,
    QuantityVector,
    cross_cartesian_vectors,
    validate_input,
    validate_output,
)


def angular_momentum_definition(
    position_vector_: Vector,
    linear_momentum_: Vector,
) -> Vector:
    r"""
    Pseudovector of angular momentum is defined as the cross product between the position vector and the momentum vector.

    Law:
        :code:`L = cross(r, p)`

    Latex:
        .. math::
            \vec L = \vec r \times \vec p

    :param position_vector\_: position vector of the particle relative to a fixed point.

        Symbol: :code:`r`

        Latex: :math:`\vec r`

        Dimension: *length*

    :param linear_momentum\_: vector of linear momentum of the particle.

        Symbol: :code:`p`

        Latex: :math:`\vec p`

        Dimension: *momentum*

    :return: pseudovector of angular momentum of the particle.

        Symbol: :code:`L`

        Latex: :math:`\vec L`

        Dimension: *length* * *momentum*
    """

    return cross_cartesian_vectors(position_vector_, linear_momentum_)


@validate_input(
    position_vector_=units.length,
    linear_momentum_=units.momentum,
)
@validate_output(units.length * units.momentum)
def calculate_angular_momentum(
    position_vector_: QuantityVector,
    linear_momentum_: QuantityVector,
) -> QuantityVector:
    result_vector = angular_momentum_definition(
        position_vector_.to_base_vector(),
        linear_momentum_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result_vector)
