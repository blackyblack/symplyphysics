from typing import Sequence
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

# Description
## The center of mass (com) of a system of particles is a unique point at any given time where
## the weighted relative position of the distributed mass sums to zero.

# Law: r_com = sum_i(m_i * r_i) / sum(m_i)
## r_com - position vector of system's center of mass
## m_i - mass of i'th particle
## r_i - position vector of i'th particle
## sum_i(x_i) denotes the sum of all x_i over index i


def center_of_mass_law(
    masses_: Sequence[ScalarValue],
    position_vectors_: Sequence[Vector],
) -> Vector:
    assert len(masses_) == len(position_vectors_)
    assert len(masses_) > 0, "At least one particle should be present."
    total_mass = 0.0
    first = True
    for mass, position_vector in zip(masses_, position_vectors_):
        total_mass += mass
        scaled = scale_vector(mass, position_vector)
        result: Vector = scaled if first else add_cartesian_vectors(result, scaled)
        first = False
    scaled_result = scale_vector(1 / total_mass, result)
    return scaled_result


@validate_input(
    masses_=units.mass,
    # position_vectors_=units.length,  # FIXME: function raises errors.UnitsError
)
@validate_output(units.length)
def calculate_center_of_mass(
    masses_: Sequence[Quantity],
    position_vectors_: Sequence[QuantityVector],
) -> QuantityVector:
    result = center_of_mass_law(masses_, position_vectors_)
    return QuantityVector(result.components, result.coordinate_system)
