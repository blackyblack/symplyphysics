from functools import reduce
from operator import add
from typing import Sequence
from sympy import S, Expr, cos, sin, sqrt, sympify
from ..coordinate_systems.coordinate_systems import CoordinateSystem
from .vectors import Vector, VectorComponent
from ..expr_comparisons import expr_equals


def equal_lists(list_left: Sequence[VectorComponent],
    list_right: Sequence[VectorComponent]) -> bool:
    max_len = max(len(list_left), len(list_right))
    list_left_extended = list(list_left) + [0] * (max_len - len(list_left))
    list_right_extended = list(list_right) + [0] * (max_len - len(list_right))
    for l, r in zip(list_left_extended, list_right_extended):
        if not expr_equals(l, r):
            return False
    return True


def add_lists(list_left: Sequence[VectorComponent],
    list_right: Sequence[VectorComponent]) -> list[Expr]:
    max_len = max(len(list_left), len(list_right))
    list_left_extended = list(list_left) + [0] * (max_len - len(list_left))
    list_right_extended = list(list_right) + [0] * (max_len - len(list_right))
    return list(map(lambda lr: sympify(lr[0] + lr[1]), zip(list_left_extended, list_right_extended)))


# Multiply each element of the list to 'scalar_value'
def scale_list(scalar_value: Expr, list_: Sequence[VectorComponent]) -> list[Expr]:
    return [scalar_value * e for e in list_]


# Multiply elements of two lists respectively and sum the results
def multiply_lists_and_sum(list_left: Sequence[VectorComponent],
    list_right: Sequence[VectorComponent]) -> Expr:
    return sympify(reduce(add, map(lambda lr: lr[0] * lr[1], zip(list_left, list_right)), 0))


# Resulting vector is orthogonal to both input vectors and its magnitude is equal to
# the area of the parallelogram spanned by input vectors
def cross_multiply_lists(list_left: Sequence[VectorComponent],
    list_right: Sequence[VectorComponent]) -> list[Expr]:
    dimensions = 3
    if len(list_left) > dimensions:
        raise ValueError(f"Cross product is only defined for {dimensions} dimensions. Got: {len(list_left)}")
    if len(list_right) > dimensions:
        raise ValueError(f"Cross product is only defined for {dimensions} dimensions. Got: {len(list_right)}")
    list_left_extended = list(list_left) + [0] * (dimensions - len(list_left))
    list_right_extended = list(list_right) + [0] * (dimensions - len(list_right))
    ax, ay, az = list_left_extended
    bx, by, bz = list_right_extended
    return [sympify(ay * bz - az * by), sympify(az * bx - ax * bz), sympify(ax * by - ay * bx)]


# Compare two Python vectors
def equal_vectors(vector_left: Vector, vector_right: Vector) -> bool:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(
            f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}"
        )
    return equal_lists(vector_left.components, vector_right.components)


# Sum of two Python vectors
# Sum of two vectors can be seen as a diagonal of the parallelogram, where vectors are adjacent sides of this parallelogram.
# To subtract vectors, multiply one of the vectors to -1 and add them.
#NOTE: adding two non-cartesian vectors is not a trivial task. We suggest to convert them to cartesian
#      vectors, add them and convert back.
def add_cartesian_vectors(vector_left: Vector, vector_right: Vector) -> Vector:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(
            f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}"
        )
    if vector_left.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(
            vector_left.coordinate_system.coord_system_type)
        raise ValueError(
            f"Addition only supported for cartesian coordinates: got {coord_name_from}")
    if vector_right.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(
            vector_right.coordinate_system.coord_system_type)
        raise ValueError(
            f"Addition only supported for cartesian coordinates: got {coord_name_from}")
    result = add_lists(vector_left.components, vector_right.components)
    return Vector(vector_left.coordinate_system, result)


# Change Vector magnitude (length)
# Scalar multiplication changes the magnitude of the vector and does not change it's direction.
def scale_vector(scalar_value: Expr, vector: Vector) -> Vector:
    if vector.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        vector_components = scale_list(scalar_value, vector.components)
        return Vector(vector.coordinate_system, vector_components)
    if vector.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        vector_components = [sympify(c) for c in vector.components]
        vector_size = len(vector_components)
        if vector_size > 0:
            vector_components[0] = vector_components[0] * scalar_value
        if vector_size > 2:
            vector_components[2] = vector_components[2] * scalar_value
        return Vector(vector.coordinate_system, vector_components)
    if vector.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        vector_components = [sympify(c) for c in vector.components]
        if len(vector_components) > 0:
            vector_components[0] = vector_components[0] * scalar_value
        return Vector(vector.coordinate_system, vector_components)
    # never
    return vector


# Dot product of two Python vectors
# Dot product equals to magnitudes of both vectors multiplied * cos(phi), where
# phi is angle between vectors.
# Hence vectors are orthogonal (perpendicular) when dot product is zero.
def dot_vectors(vector_left: Vector, vector_right: Vector) -> Expr:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(
            f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}"
        )
    if vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        return multiply_lists_and_sum(vector_left.components, vector_right.components)
    if vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        left_components = list(vector_left.components)
        left_components.extend([0] * (3 - len(left_components)))
        right_components = list(vector_right.components)
        right_components.extend([0] * (3 - len(right_components)))
        r1, r2 = left_components[0], right_components[0]
        theta1, theta2 = left_components[1], right_components[1]
        z1, z2 = left_components[2], right_components[2]
        return r1 * r2 * cos(theta1 - theta2) + z1 * z2
    if vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        left_components = list(vector_left.components)
        left_components.extend([0] * (3 - len(left_components)))
        right_components = list(vector_right.components)
        right_components.extend([0] * (3 - len(right_components)))
        r1, r2 = left_components[0], right_components[0]
        theta1, theta2 = left_components[1], right_components[1]
        phi1, phi2 = left_components[2], right_components[2]
        return r1 * r2 * (sin(theta1) * sin(theta2) * cos(phi1 - phi2) + cos(theta1) * cos(theta2))
    # never
    return S.Zero


def vector_magnitude(vector_: Vector | Sequence[VectorComponent]) -> Expr:
    squared_sum = dot_vectors(vector_, vector_) if isinstance(vector_,
        Vector) else multiply_lists_and_sum(vector_, vector_)
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
            f"Cross product only supported for cartesian coordinates: got {coord_name_from}")
    if vector_right.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(
            vector_right.coordinate_system.coord_system_type)
        raise ValueError(
            f"Cross product only supported for cartesian coordinates: got {coord_name_from}")
    result = cross_multiply_lists(vector_left.components, vector_right.components)
    return Vector(vector_left.coordinate_system, result)
