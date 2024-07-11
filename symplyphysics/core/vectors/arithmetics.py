from functools import reduce
from operator import add
from typing import Optional, Sequence
from sympy import S, Expr, cos, sin, sqrt, sympify, diff, integrate

from .vectors import Vector
from ..expr_comparisons import expr_equals
from ..coordinate_systems.coordinate_systems import CoordinateSystem
from ..dimensions import ScalarValue


# Add zeroes so that both vectors have the same length.
# Use 'max_size' to increase or trim vector size. Vectors are aligned to the
# larger (more dimensions) vector size if not set.
def _extend_two_vectors(
        vector_left: Vector,
        vector_right: Vector,
        max_size: Optional[int] = None) -> tuple[Sequence[ScalarValue], Sequence[ScalarValue]]:
    max_size = max(len(vector_left.components), len(
        vector_right.components)) if max_size is None else max_size
    list_left_extended = list(
        vector_left.components) + [0] * (max_size - len(vector_left.components))
    list_right_extended = list(
        vector_right.components) + [0] * (max_size - len(vector_right.components))
    return (list_left_extended, list_right_extended)


# Compare two vectors
def equal_vectors(vector_left: Vector, vector_right: Vector) -> bool:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(
            f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}"
        )
    (list_left_extended, list_right_extended) = _extend_two_vectors(vector_left, vector_right)
    for l, r in zip(list_left_extended, list_right_extended):
        if not expr_equals(l, r):
            return False
    return True


# Add any number of vectors
#NOTE: adding two non-cartesian vectors is not a trivial task. We suggest to convert them to cartesian
#      vectors, add them and convert back.
def add_cartesian_vectors(*vectors: Vector) -> Vector:
    def add_two_cartesian_vectors(vector_left: Vector, vector_right: Vector) -> Vector:
        if vector_left.coordinate_system != vector_right.coordinate_system:
            raise ValueError(
                "Vectors must have the same coordinate system, "
                f"got {vector_left.coordinate_system} and {vector_right.coordinate_system}"
            )

        left_type = vector_left.coordinate_system.coord_system_type
        right_type = vector_right.coordinate_system.coord_system_type

        if left_type != CoordinateSystem.System.CARTESIAN or right_type != CoordinateSystem.System.CARTESIAN:
            left_name = CoordinateSystem.system_to_transformation_name(left_type)
            right_name = CoordinateSystem.system_to_transformation_name(right_type)
            raise ValueError(f"Expect two Cartesian vectors, got {left_name} and {right_name}")

        (list_left_extended, list_right_extended) = _extend_two_vectors(vector_left, vector_right)
        result = [
            sympify(left_component + right_component)
            for (left_component, right_component)
            in zip(list_left_extended, list_right_extended)
        ]

        return Vector(result, vector_left.coordinate_system)

    if len(vectors) < 1:
        raise ValueError("Provide at least one vector")

    result = vectors[0]
    for vector in vectors[1:]:
        result = add_two_cartesian_vectors(result, vector)
    return result


# Subtract vectors from the first given vector
def subtract_cartesian_vectors(*vectors: Vector) -> Vector:
    if len(vectors) < 2:
        raise ValueError("Provide at least two vectors")

    vector_minuend = vectors[0]
    vector_subtrahend = add_cartesian_vectors(*vectors[1:])
    return add_cartesian_vectors(
        vector_minuend,
        scale_vector(-1, vector_subtrahend),
    )


# Change Vector magnitude (length)
# Scalar multiplication changes the magnitude of the vector and does not change it's direction.
def scale_vector(scalar_value: ScalarValue, vector: Vector) -> Vector:
    if vector.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        vector_components = [scalar_value * e for e in vector.components]
        return Vector(vector_components, vector.coordinate_system)
    if vector.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        vector_components = list(vector.components)
        vector_size = len(vector_components)
        if vector_size > 0:
            vector_components[0] = vector_components[0] * scalar_value
        if vector_size > 2:
            vector_components[2] = vector_components[2] * scalar_value
        return Vector(vector_components, vector.coordinate_system)
    if vector.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        vector_components = list(vector.components)
        if len(vector_components) > 0:
            vector_components[0] = vector_components[0] * scalar_value
        return Vector(vector_components, vector.coordinate_system)
    # never
    return vector


# Multiply elements of two lists respectively and sum the results
def _multiply_lists_and_sum(list_left: Sequence[ScalarValue],
    list_right: Sequence[ScalarValue]) -> Expr:
    return sympify(reduce(add, map(lambda lr: lr[0] * lr[1], zip(list_left, list_right)), 0))


# Dot product of two vectors
# Dot product equals to magnitudes of both vectors multiplied * cos(phi), where
# phi is angle between vectors.
# Hence vectors are orthogonal (perpendicular) when dot product is zero.
def dot_vectors(vector_left: Vector, vector_right: Vector) -> Expr:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(
            f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}"
        )
    dimensions = 3
    if vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        return _multiply_lists_and_sum(vector_left.components, vector_right.components)
    if vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        (list_left_extended,
            list_right_extended) = _extend_two_vectors(vector_left, vector_right, dimensions)
        r1, theta1, z1 = list_left_extended
        r2, theta2, z2 = list_right_extended
        return r1 * r2 * cos(theta1 - theta2) + z1 * z2
    if vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        (list_left_extended,
            list_right_extended) = _extend_two_vectors(vector_left, vector_right, dimensions)
        r1, theta1, phi1 = list_left_extended
        r2, theta2, phi2 = list_right_extended
        return r1 * r2 * (sin(phi1) * sin(phi2) * cos(theta1 - theta2) + cos(phi1) * cos(phi2))
    # never
    return S.Zero


def vector_magnitude(vector_: Vector) -> Expr:
    squared_sum = dot_vectors(vector_, vector_)
    return sqrt(squared_sum)


#NOTE: vector multiplication of two non-cartesian vectors is not a trivial task. We suggest to convert them to cartesian
#      vectors, multiply them and convert back.
#NOTE: cross product with all 3 basic properties:
#      - Cross product is bilinear function of input vectors;
#      - Resulting vector is perpendicular to both input vectors;
#      - Magnitude of resulting vector equals to the area of the parallelogram spanned by input vectors;
#      is defined only in 3 and 7 dimensional Euclidean space.
def cross_cartesian_vectors(vector_left: Vector, vector_right: Vector) -> Vector:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(
            f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}"
        )
    if vector_left.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(
            vector_left.coordinate_system.coord_system_type)
        raise ValueError(
            f"Cross product is only supported for cartesian coordinates: got {coord_name_from}")
    if vector_right.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(
            vector_right.coordinate_system.coord_system_type)
        raise ValueError(
            f"Cross product is only supported for cartesian coordinates: got {coord_name_from}")
    dimensions = 3
    if len(vector_left.components) > dimensions:
        raise ValueError(
            f"Cross product is only defined for {dimensions} dimensions. Got: {len(vector_left.components)}"
        )
    if len(vector_right.components) > dimensions:
        raise ValueError(
            f"Cross product is only defined for {dimensions} dimensions. Got: {len(vector_right.components)}"
        )
    (list_left_extended, list_right_extended) = _extend_two_vectors(vector_left, vector_right,
        dimensions)
    ax, ay, az = list_left_extended
    bx, by, bz = list_right_extended
    result = [ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx]
    return Vector(result, vector_left.coordinate_system)


# Make unit vector (vector of size 1 and same direction as original vector)
def vector_unit(vector_: Vector) -> Vector:
    return scale_vector(1 / vector_magnitude(vector_), vector_)


# Project `original_vector_` onto `target_vector_`. The result is the orthogonal projection of
# `original_vector_` onto a straight line parallel to `target_vector_` and is parallel to 
# `target_vector_`.
def project_vector(
    original_vector_: Vector,
    target_vector_: Vector,
) -> Vector:
    return scale_vector(
        dot_vectors(original_vector_, target_vector_) / dot_vectors(target_vector_, target_vector_),
        target_vector_,
    )


# Reject `original_vector_` from `target_vector_`. The result is the orthogonal projection of
# `original_vector_` onto the (hyper)plane orthogonal to `target_vector_` and is orthogonal to
# `target_vector_`.
def reject_cartesian_vector(
    original_vector_: Vector,
    target_vector_: Vector,
) -> Vector:
    return add_cartesian_vectors(
        original_vector_,
        scale_vector(-1, project_vector(original_vector_, target_vector_))
    )


# Only works in non-rotating coordinate systems
def diff_cartesian_vector(
    vector_: Vector,
    *args: Expr,
) -> Vector:
    if vector_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Component-wise vector differentiation is only supported for Cartesian coordinates")

    components = [diff(component, *args) for component in vector_.components]

    return Vector(components, vector_.coordinate_system)


def integrate_cartesian_vector(
    vector_: Vector,
    *args: Expr | tuple[Expr, Expr, Expr],
) -> Vector:
    if vector_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Component-wise vector integration is only supported for Cartesian coordinates")

    components = [integrate(component, *args) for component in vector_.components]

    return Vector(components, vector_.coordinate_system)
