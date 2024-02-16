from typing import Sequence

from sympy import S
from symplyphysics import (
    CoordinateSystem,
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
    for c in position_vectors_:
        if c.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
            raise ValueError(
                f"Radius vector {c} should be in cartesian coordinate system"
            )

    position_base_vectors = [v.to_base_vector() for v in position_vectors_]
    result_vector = center_of_mass_law(masses_, position_base_vectors)
    return QuantityVector.from_base_vector(result_vector)
