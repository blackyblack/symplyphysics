from typing import Sequence
from symplyphysics import (units, validate_output)
from symplyphysics.core.dimensions import assert_equivalent_dimension
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
    forces_ = list(forces_)
    result = Vector([]) if len(forces_) == 0 else Vector([], forces_[0].coordinate_system)
    for f in forces_:
        result = add_cartesian_vectors(result, f)
    return result


@validate_output(units.force)
def calculate_resultant_force(forces_: Sequence[QuantityVector]) -> QuantityVector:
    if len(forces_) == 0:
        return QuantityVector([])
    forces_ = list(forces_)
    for f in forces_:
        assert_equivalent_dimension(f.dimension, f.display_name, "calculate_resultant_force",
            units.force)
    return QuantityVector(superposition_law(forces_).components, forces_[0].coordinate_system)
