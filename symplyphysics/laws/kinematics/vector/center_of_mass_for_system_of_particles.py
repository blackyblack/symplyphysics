"""
Center of mass for a system of particles
========================================

The center of mass (com) of a system of particles is a unique point at any given time where
the sum of weighted relative positions of the distributed mass is zero.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Center_of_mass#A_system_of_particles>`__.
"""

from typing import Sequence
from sympy import S
from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Quantity,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
    scale_vector,
)
from symplyphysics.core.dimensions import ScalarValue


def center_of_mass_law(
    masses_: Sequence[ScalarValue],
    position_vectors_: Sequence[Vector],
) -> Vector:
    r"""
    Vector of the center of mass from masses and position vectors.

    Law:
        :code:`r_com = Sum(m_i * r_i, i) / Sum(m_i, i)`

    Latex:
        .. math::
            {\vec r}_\text{com} = \frac{\sum_i m_i {\vec r}_i}{\sum_i m_i}

    :param masses\_: sequence of masses of individual parts

        Symbol: :code:`m_i`

        Latex: :math:`m_i`

        Dimension: *mass*

    :param position_vectors\_: sequence of position vectors of individual parts

        Symbol: :code:`r_i`

        Latex: :math:`{\vec r}_i`

        Dimension: *length*

    :return: vector of the center of mass

        Symbol: :code:`r_com`

        Latex: :math:`{\vec r}_\text{com}`

        Dimension: *length*
    """

    result = Vector([0, 0, 0], next(iter(position_vectors_)).coordinate_system)
    total_mass = S.Zero
    for mass, position_vector in zip(masses_, position_vectors_):
        total_mass += mass
        scaled = scale_vector(mass, position_vector)
        result = add_cartesian_vectors(result, scaled)
    scaled_result = scale_vector(1 / total_mass, result)
    return scaled_result


@validate_input(
    masses_=units.mass,
    position_vectors_=units.length,
)
@validate_output(units.length)
def calculate_center_of_mass(
    masses_: Sequence[Quantity],
    position_vectors_: Sequence[QuantityVector],
) -> QuantityVector:
    if len(masses_) != len(position_vectors_):
        raise ValueError("Mass and position arrays should have the same lengths")
    if len(position_vectors_) == 0:
        raise ValueError("At least one particle should be present")
    position_base_vectors = [v.to_base_vector() for v in position_vectors_]
    result_vector = center_of_mass_law(masses_, position_base_vectors)
    return QuantityVector.from_base_vector(result_vector)
