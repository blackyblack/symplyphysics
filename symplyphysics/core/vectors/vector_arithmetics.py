from typing import List
from sympy import Expr
from ..expr_comparisons import expr_equals


# Compare two Python vectors
def equal_vectors(vector_left: List[Expr], vector_right: List[Expr]) -> bool:
    for i in range(max(len(vector_left), len(vector_right))):
        val1 = 0 if i >= len(vector_left) else vector_left[i]
        val2 = 0 if i >= len(vector_right) else vector_right[i]
        if not expr_equals(val1, val2): return False
    return True

# Sum of two Python vectors
def add_vectors(vector_left: List[Expr], vector_right: List[Expr]) -> List[Expr]:
    result = []
    for i in range(max(len(vector_left), len(vector_right))):
        val1 = 0 if i >= len(vector_left) else vector_left[i]
        val2 = 0 if i >= len(vector_right) else vector_right[i]
        result.append(val1 + val2)
    return result

# Scalar Python vector multiplication
def multiply_vector(scalar_value: Expr, vector: List[Expr]) -> List[Expr]:
    return [scalar_value * e for e in vector]

# Dot product of two Python vectors
def dot_vectors(vector_left: List[Expr], vector_right: List[Expr]) -> Expr:
    result = 0
    for i in range(min(len(vector_left), len(vector_right))):
        result += (vector_left[i] * vector_right[i])
    return result
