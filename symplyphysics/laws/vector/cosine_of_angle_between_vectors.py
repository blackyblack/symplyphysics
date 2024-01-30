from sympy import Eq, Symbol as SymSymbol, cos, solve
from symplyphysics import (
    Quantity,
    print_expression,
    Vector,
    dot_vectors,
    vector_magnitude,
)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.symbols.quantities import scale_factor

# Description:
## The angle between two vectors (in any dimensions) can be found via their dot product
## and their lengths.

# Law: cos(phi) = (a, b) / (norm(a) * norm(b))
## a, b - vectors
## phi - angle between vectors a and b
## norm(v) - norm, or length, of vector v

angle_between_vectors = SymSymbol("cosine_between_vectors", real=True)
dot_product = SymSymbol("dot_product", real=True)
vector_left_norm = SymSymbol("vector_left_norm", positive=True)
vector_right_norm = SymSymbol("vector_right_norm", positive=True)

law = Eq(cos(angle_between_vectors), dot_product / (vector_left_norm * vector_right_norm))


def print_law() -> str:
    return print_expression(law)


def cosine_law(vector_left_: Vector, vector_right_: Vector) -> ScalarValue:
    dot_product_ = dot_vectors(vector_left_, vector_right_)
    vector_left_norm_ = vector_magnitude(vector_left_)
    vector_right_norm_ = vector_magnitude(vector_right_)
    return law.rhs.subs({
        dot_product: dot_product_,
        vector_left_norm: vector_left_norm_,
        vector_right_norm: vector_right_norm_,
    })


def calculate_angle_between_vectors(
    dot_product_: Quantity | float,
    vector_left_norm_: Quantity | float,
    vector_right_norm_: Quantity | float,
) -> float:
    dot_product_ = scale_factor(dot_product_)
    vector_left_norm_ = scale_factor(vector_left_norm_)
    vector_right_norm_ = scale_factor(vector_right_norm_)
    result = solve(law, angle_between_vectors)[1].subs({
        dot_product: dot_product_,
        vector_left_norm: vector_left_norm_,
        vector_right_norm: vector_right_norm_,
    })
    return result
