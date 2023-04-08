from sympy import Expr, cos, sin, sqrt
from ..coordinate_systems.coordinate_systems import CoordinateSystem
from .vectors import Vector
from ..expr_comparisons import expr_equals


# Compare two Python vectors
def equal_vectors(vector_left: Vector, vector_right: Vector) -> bool:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}")
    for i in range(max(len(vector_left.components), len(vector_right.components))):
        val1 = 0 if i >= len(vector_left.components) else vector_left.components[i]
        val2 = 0 if i >= len(vector_right.components) else vector_right.components[i]
        if not expr_equals(val1, val2): return False
    return True

# Sum of two Python vectors
# Sum of two vectors can be seen as a diagonal of the parallelogram, where vectors are adjacent sides of this parallelogram.
# To subtract vectors, multiply one of the vectors to -1 and add them.
#NOTE: adding two non-cartesian vectors is not a trivial task. We suggest to convert them to cartesian
#      vectors, add them and convert back.
def add_cartesian_vectors(vector_left: Vector, vector_right: Vector) -> Vector:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}")
    if vector_left.coordinate_system is not None and vector_left.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(vector_left.coordinate_system.coord_system_type)
        raise ValueError(f"Addition only supported for cartesian coordinates: got {coord_name_from}")
    if vector_right.coordinate_system is not None and vector_right.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(vector_right.coordinate_system.coord_system_type)
        raise ValueError(f"Addition only supported for cartesian coordinates: got {coord_name_from}")
    result = []
    for i in range(max(len(vector_left.components), len(vector_right.components))):
        val1 = 0 if i >= len(vector_left.components) else vector_left.components[i]
        val2 = 0 if i >= len(vector_right.components) else vector_right.components[i]
        result.append(val1 + val2)
    return Vector(result, vector_left.coordinate_system)

# Change Vector magnitude (length)
# Scalar multiplication changes the magnitude of the vector and does not change it's direction.
def scale_vector(scalar_value: Expr, vector: Vector) -> Vector:
    vector_size = len(vector.components)
    if vector.coordinate_system is None or vector.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        vector_components = [scalar_value * e for e in vector.components]
        return Vector(vector_components, vector.coordinate_system)
    if vector.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        vector_components = [*vector.components]
        if vector_size > 0:
            vector_components[0] = vector_components[0] * scalar_value
        if vector_size > 2:
            vector_components[2] = vector_components[2] * scalar_value
        return Vector(vector_components, vector.coordinate_system)
    if vector.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        vector_components = [*vector.components]
        if vector_size > 0:
            vector_components[0] = vector_components[0] * scalar_value
        return Vector(vector_components, vector.coordinate_system)
    coord_name_from = CoordinateSystem.system_to_transformation_name(vector.coordinate_system.coord_system_type)
    raise ValueError(f"Scaling only supported for cartesian, cylindrical or spherical coordinates: got {coord_name_from}")

# Dot product of two Python vectors
# Dot product equals to magnitudes of both vectors multiplied * cos(phi), where
# phi is angle between vectors.
# Hence vectors are orthogonal (perpendicular) when dot product is zero.
def dot_vectors(vector_left: Vector, vector_right: Vector) -> Expr:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}")
    
    if vector_left.coordinate_system is None or vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        result = 0
        for i in range(min(len(vector_left.components), len(vector_right.components))):
            result += (vector_left.components[i] * vector_right.components[i])
        return result
    if vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        left_components = vector_left.components
        left_components.extend([0] * (3 - len(left_components)))
        right_components = vector_right.components
        right_components.extend([0] * (3 - len(right_components)))
        r1, r2 = left_components[0], right_components[0]
        theta1, theta2 = left_components[1], right_components[1]
        z1, z2 = left_components[2], right_components[2]
        return r1 * r2 * cos(theta1 - theta2) + z1 * z2
                       
    if vector_left.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        left_components = vector_left.components
        left_components.extend([0] * (3 - len(left_components)))
        right_components = vector_right.components
        right_components.extend([0] * (3 - len(right_components)))
        r1, r2 = left_components[0], right_components[0]
        theta1, theta2 = left_components[1], right_components[1]
        phi1, phi2 = left_components[2], right_components[2]
        x = r1 * r2 * (sin(theta1) * sin(theta2) * cos(phi1 - phi2) + cos(theta1) * cos(theta2))
        return x
    
def vector_magnitude(vector_: Vector) -> Expr:
    return sqrt(dot_vectors(vector_, vector_))

#NOTE: vector multiplication of two non-cartesian vectors is not a trivial task. We suggest to convert them to cartesian
#      vectors, multiply them and convert back.
#NOTE: cross product with all 3 basic properties:
#      - Cross product is bilinear function of input vectors;
#      - Resulting vector is perpendicular to both input vectors;
#      - Magnitude of resulting vector equals to the area of the parallelogram spanned by input vectors;
#      is defined only in 3 and 7 dimensional Euclidean space.
def cross_cartesian_vectors(vector_left: Vector, vector_right: Vector) -> Vector:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}")
    if vector_left.coordinate_system is not None and vector_left.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(vector_left.coordinate_system.coord_system_type)
        raise ValueError(f"Cross product only supported for cartesian coordinates: got {coord_name_from}")
    if vector_right.coordinate_system is not None and vector_right.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(vector_right.coordinate_system.coord_system_type)
        raise ValueError(f"Cross product only supported for cartesian coordinates: got {coord_name_from}")
    if len(vector_left.components) > 3:
        raise ValueError(f"Cross product is only defined for 3 dimensions. Got: {len(vector_left.components)}")
    if len(vector_right.components) > 3:
        raise ValueError(f"Cross product is only defined for 3 dimensions. Got: {len(vector_right.components)}")
    
    left_components = vector_left.components
    left_components.extend([0] * (3 - len(left_components)))
    right_components = vector_right.components
    right_components.extend([0] * (3 - len(right_components)))
    ax, ay, az = left_components
    bx, by, bz = right_components
    return Vector([ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx], vector_left.coordinate_system)
