from sympy import Eq, Symbol as SymSymbol, sin, sqrt
from symplyphysics import (
    print_expression,
    Vector,
    vector_magnitude,
)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.laws.vector import cosine_of_angle_between_vectors as cosine_law

# Description
## The norm of the cross product vector between two vectors is the product of their
## respective norms and the sine of the angle between them.

# Law: norm(c) = norm(a) * norm(b) * sin(phi)
## a, b - vectors in 3D
## phi - angle between vectors a and b
## c = [a, b] - cross product of vectors a and b
## norm(v) - norm, or length, of vector v

# Notes
## - The direction of the cross product vector can be found via the right hand rule

cross_product_norm = SymSymbol("cross_product_norm", positive=True)
vector_left_norm = SymSymbol("vector_left_norm", positive=True)
vector_right_norm = SymSymbol("vector_right_norm", positive=True)
angle_between_vectors = SymSymbol("angle_between_vectors", real=True)

law = Eq(cross_product_norm, vector_left_norm * vector_right_norm * sin(angle_between_vectors))


def print_law() -> str:
    return print_expression(law)


def cross_product_norm_law(vector_left_: Vector, vector_right_: Vector) -> ScalarValue:
    if any(len(v.components) != 3 for v in (vector_left_, vector_right_)):
        raise ValueError("Cross product is calculated only for 3-dimensional spaces.")

    vector_left_norm_ = vector_magnitude(vector_left_)
    vector_right_norm_ = vector_magnitude(vector_right_)
    cos_between_vectors = cosine_law.cosine_law(vector_left_, vector_right_)
    sin_between_vectors = sqrt(1 - cos_between_vectors**2)
    return law.rhs.subs({
        vector_left_norm: vector_left_norm_,
        vector_right_norm: vector_right_norm_,
        sin(angle_between_vectors): sin_between_vectors,
    })
