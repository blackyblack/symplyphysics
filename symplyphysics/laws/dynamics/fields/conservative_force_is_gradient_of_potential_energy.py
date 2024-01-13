from symplyphysics import Vector
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.fields.operators import gradient_operator

# Description
## A conservative force is a such a force, the total work of which in moving a particle
## between two points is independent of the path taken. Alternative definition states that
## if a particle travels in a closed loop, the total work done by a conservative force is zero.

## For every conservative force there exists a scalar field (potential), for which applies
## the following equality:

# Law: F = - grad U
## F - vector field of the conservative force
## U - scalar field of its potential
## grad - gradient operator

def law(potential_: ScalarField) -> VectorField:
    gradient_field = gradient_operator(potential_).apply_to_basis()
    result_vector = Vector(
        [-component for component in gradient_field.components],
        gradient_field.coordinate_system
    )
    return VectorField.from_vector(result_vector)
