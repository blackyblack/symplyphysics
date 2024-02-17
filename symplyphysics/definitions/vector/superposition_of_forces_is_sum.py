from typing import Sequence
from symplyphysics import (CoordinateSystem, units, validate_input, validate_output)
from symplyphysics.core.vectors.arithmetics import add_cartesian_vectors
from symplyphysics.core.vectors.vectors import QuantityVector, Vector

# Description
## R = sum(F)
## Where:
## F is one of the force vectors, acting on the body,
## R - resultant (or net) force vector.

# Conditions:
## - Vectors are defined in terms of cartesian coordinate system.


def superposition_law(forces_: Sequence[Vector]) -> Vector:
    result = Vector([0, 0, 0], next(iter(forces_)).coordinate_system)
    for f in forces_:
        result = add_cartesian_vectors(result, f)
    return result


@validate_input(forces_=units.force)
@validate_output(units.force)
def calculate_resultant_force(forces_: Sequence[QuantityVector]) -> QuantityVector:
    if len(forces_) == 0:
        raise ValueError("At least one force should be present")
    for c in forces_:
        if c.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
            raise ValueError(f"Force vector {c} should be in cartesian coordinate system")
    forces_base_vectors = [f.to_base_vector() for f in forces_]
    result_force_vector = superposition_law(forces_base_vectors)
    return QuantityVector.from_base_vector(result_force_vector)
