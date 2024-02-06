from sympy import Eq, Symbol as SymSymbol, cos, solve
from symplyphysics import Quantity, print_expression
from symplyphysics.core.symbols.quantities import scale_factor

# Description:
## Dot product is scalar binary operation defined as the product of the norm of the vectors
## and the cosine of the angle between them.

# Law: dot(a, b) = norm(a) * norm(b) * cos(phi)
## a, b - vectors
## phi - angle between vectors a and b
## dot(a, b) - dot product between vectors a and b
## norm(v) - norm, or length, of vector v

dot_product = SymSymbol("dot_product", real=True)
vector_left_norm = SymSymbol("vector_left_norm", positive=True)
vector_right_norm = SymSymbol("vector_right_norm", positive=True)
angle_between_vectors = SymSymbol("cosine_between_vectors", real=True)

law = Eq(dot_product, vector_left_norm * vector_right_norm * cos(angle_between_vectors))


def print_law() -> str:
    return print_expression(law)


def calculate_cosine_between_vectors(
    dot_product_: Quantity | float,
    vector_left_norm_: Quantity | float,
    vector_right_norm_: Quantity | float,
) -> Quantity:
    dot_product_ = scale_factor(dot_product_)
    vector_left_norm_ = scale_factor(vector_left_norm_)
    vector_right_norm_ = scale_factor(vector_right_norm_)
    result = solve(law, cos(angle_between_vectors))[0].subs({
        dot_product: dot_product_,
        vector_left_norm: vector_left_norm_,
        vector_right_norm: vector_right_norm_,
    })
    return Quantity(result)
