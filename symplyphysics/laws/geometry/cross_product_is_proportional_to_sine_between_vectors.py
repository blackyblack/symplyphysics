from sympy import Eq, Symbol as SymSymbol, sin, solve
from symplyphysics import Quantity
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## Cross product is a binary vector operation, defined as such:
## - Its magnitude is the product of the magnitudes of the two vectors and the sine between them.
## - Its direction can be found via the right-hand rule.

# Law: norm(cross(a, b)) = norm(a) * norm(b) * sin(phi)
## a, b - vectors
## phi - angle between them
## cross(a, b) - cross product between a and b
## norm(v) - norm, or length, of a

# Condition: cross product is only defined for 3- (and 7-) dimensional vectors.

# TODO: update documentation

cross_product_norm = SymSymbol("cross_product_norm", positive=True)
vector_left_norm = SymSymbol("vector_left_norm", positive=True)
vector_right_norm = SymSymbol("vector_right_norm", positive=True)
angle_between_vectors = SymSymbol("cosine_between_vectors", real=True)

law = Eq(cross_product_norm, vector_left_norm * vector_right_norm * sin(angle_between_vectors))


def calculate_sine_between_vectors(
    cross_product_norm_: Quantity | float,
    vector_left_norm_: Quantity | float,
    vector_right_norm_: Quantity | float,
) -> Quantity:
    cross_product_norm_ = scale_factor(cross_product_norm_)
    vector_left_norm_ = scale_factor(vector_left_norm_)
    vector_right_norm_ = scale_factor(vector_right_norm_)
    result = solve(law, sin(angle_between_vectors))[0].subs({
        cross_product_norm: cross_product_norm_,
        vector_left_norm: vector_left_norm_,
        vector_right_norm: vector_right_norm_,
    })
    return Quantity(result)
