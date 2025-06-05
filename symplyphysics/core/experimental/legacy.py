from sympy import Expr

from ..vectors.vectors import Vector as LegacyVector
from ..vectors.arithmetics import (add_cartesian_vectors as legacy_add_cartesian_vectors,
    scale_vector as legacy_scale_vector)

from .coordinate_systems import CoordinateVector, CartesianCoordinateSystem, CARTESIAN
from .vectors import into_terms, split_factor


def into_legacy_vector(expr: Expr) -> LegacyVector:
    result = LegacyVector([0, 0, 0])

    for term in into_terms(expr):
        vector, factor = split_factor(term)

        if not isinstance(vector, CoordinateVector):
            raise TypeError(f"Expected CoordinateVector, got {type(term).__name__}")

        if not isinstance(vector.system, CartesianCoordinateSystem):
            raise TypeError(f"Unsupported system type {type(vector.system).__name__}")

        legacy_vector = LegacyVector(vector.components)

        legacy_term = legacy_scale_vector(factor, legacy_vector)

        result = legacy_add_cartesian_vectors(result, legacy_term)

    return result


def from_legacy_vector(vector: LegacyVector) -> Expr:
    return CoordinateVector(vector.components, CARTESIAN)
