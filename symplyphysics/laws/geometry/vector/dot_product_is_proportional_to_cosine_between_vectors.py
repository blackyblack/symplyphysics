from sympy import solve, cos
from symplyphysics import (
    Quantity,
    QuantityVector,
    Vector,
    dot_vectors,
    vector_magnitude,
)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.laws.geometry import (
    dot_product_is_proportional_to_cosine_between_vectors as cosine_scalar_law,
)

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
    return solve(
        cosine_scalar_law.law, cos(cosine_scalar_law.angle_between_vectors)
    )[0].subs({
        cosine_scalar_law.dot_product: dot_product_,
        cosine_scalar_law.vector_left_norm: vector_left_norm_,
        cosine_scalar_law.vector_right_norm: vector_right_norm_,
    })


def calculate_cosine_between_vectors(
    vector_left_: QuantityVector, vector_right_: QuantityVector
) -> Quantity:
    return Quantity(cosine_between_vectors_law(vector_left_, vector_right_))
