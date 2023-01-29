from sympy import Expr
from .vectors import Vector
from ..expr_comparisons import expr_equals


## Only works for Cartesian coordinates!
#TODO: rework to support other coordinate systems.
#TODO: add vector_length implementation

# Compare two Python vectors
def equal_vectors(vector_left: Vector, vector_right: Vector) -> bool:
    if vector_left.coord_system != vector_right.coord_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coord_system)} vs {str(vector_right.coord_system)}")
    for i in range(max(len(vector_left.components), len(vector_right.components))):
        val1 = 0 if i >= len(vector_left.components) else vector_left.components[i]
        val2 = 0 if i >= len(vector_right.components) else vector_right.components[i]
        if not expr_equals(val1, val2): return False
    return True

# Sum of two Python vectors
# Sum of two vectors can be seen as a diagonal of the parallelogram, where vectors are adjacent sides of this parallelogram.
# To subtract vectors, multiply one of the vectors to -1 and add them.
def add_vectors(vector_left: Vector, vector_right: Vector) -> Vector:
    if vector_left.coord_system != vector_right.coord_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coord_system)} vs {str(vector_right.coord_system)}")
    result = []
    for i in range(max(len(vector_left.components), len(vector_right.components))):
        val1 = 0 if i >= len(vector_left.components) else vector_left.components[i]
        val2 = 0 if i >= len(vector_right.components) else vector_right.components[i]
        result.append(val1 + val2)
    return Vector(result, vector_left.coord_system)

# Scalar Python vector multiplication
# Scalar multiplication changes the magnitude of the vector and does not change it's direction.
def multiply_vector(scalar_value: Expr, vector: Vector) -> Vector:
    vector_components = [scalar_value * e for e in vector.components]
    return Vector(vector_components, vector.coord_system)

# Dot product of two Python vectors
# Dot product equals to magnitudes of both vectors multiplied * cos(phi), where
# phi is angle between vectors.
# Hence vectors are orthogonal (perpendicular) when dot product is zero.
def dot_vectors(vector_left: Vector, vector_right: Vector) -> Expr:
    if vector_left.coord_system != vector_right.coord_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coord_system)} vs {str(vector_right.coord_system)}")
    result = 0
    for i in range(min(len(vector_left.components), len(vector_right.components))):
        result += (vector_left.components[i] * vector_right.components[i])
    return result
