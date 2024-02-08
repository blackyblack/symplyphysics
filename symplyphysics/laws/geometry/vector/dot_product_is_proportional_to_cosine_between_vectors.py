from symplyphysics import (
    Quantity,
    QuantityVector,
    Vector,
    dot_vectors,
    vector_magnitude,
)
from symplyphysics.core.dimensions import ScalarValue

# Description:
## Dot product is scalar binary operation defined as the product of the norm of the vectors
## and the cosine of the angle between them.

# Law: dot(a, b) = norm(a) * norm(b) * cos(phi)
## a, b - vectors
## phi - angle between vectors a and b
## dot(a, b) - dot product between vectors a and b
## norm(v) - norm, or length, of vector v


def cosine_between_vectors_law(
    vector_left_: Vector, vector_right_: Vector
) -> ScalarValue:
    dot_product_ = dot_vectors(vector_left_, vector_right_)
    vector_left_norm_ = vector_magnitude(vector_left_)
    vector_right_norm_ = vector_magnitude(vector_right_)
    return dot_product_ / (vector_left_norm_ * vector_right_norm_)


def calculate_cosine_between_vectors(
    vector_left_: QuantityVector, vector_right_: QuantityVector
) -> Quantity:
    return Quantity(cosine_between_vectors_law(vector_left_, vector_right_))
