from sympy import Expr

from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from .vectors import Vector
from ..expr_comparisons import expr_equals


## Only works for Cartesian coordinates!
#TODO: rework to support other coordinate systems.
#TODO: add vector_magnitude implementation
#TODO: add cross_product implementation

# Compare two Python vectors
def equal_vectors(vector_left: Vector, vector_right: Vector) -> bool:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}")
    for i in range(max(len(vector_left.components), len(vector_right.components))):
        val1 = 0 if i >= len(vector_left.components) else vector_left.components[i]
        val2 = 0 if i >= len(vector_right.components) else vector_right.components[i]
        if not expr_equals(val1, val2): return False
    return True

#TODO: does not work in polar coordinates
# Sum of two Python vectors
# Sum of two vectors can be seen as a diagonal of the parallelogram, where vectors are adjacent sides of this parallelogram.
# To subtract vectors, multiply one of the vectors to -1 and add them.
def add_vectors(vector_left: Vector, vector_right: Vector) -> Vector:
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

#TODO: does not work in polar coordinates
# Change Vector magnitude (length)
# Scalar multiplication changes the magnitude of the vector and does not change it's direction.
def scale_vector(scalar_value: Expr, vector: Vector) -> Vector:
    if vector.coordinate_system is not None and vector.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(vector.coordinate_system.coord_system_type)
        raise ValueError(f"Scaling only supported for cartesian coordinates: got {coord_name_from}")
    
    vector_components = [scalar_value * e for e in vector.components]
    return Vector(vector_components, vector.coordinate_system)

# Dot product of two Python vectors
# Dot product equals to magnitudes of both vectors multiplied * cos(phi), where
# phi is angle between vectors.
# Hence vectors are orthogonal (perpendicular) when dot product is zero.
def dot_vectors(vector_left: Vector, vector_right: Vector) -> Expr:
    if vector_left.coordinate_system != vector_right.coordinate_system:
        raise TypeError(f"Different coordinate systems in vectors: {str(vector_left.coordinate_system)} vs {str(vector_right.coordinate_system)}")
    if vector_left.coordinate_system is not None and vector_left.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(vector_left.coordinate_system.coord_system_type)
        raise ValueError(f"Dot product only supported for cartesian coordinates: got {coord_name_from}")
    if vector_right.coordinate_system is not None and vector_right.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(vector_right.coordinate_system.coord_system_type)
        raise ValueError(f"Dot product only supported for cartesian coordinates: got {coord_name_from}")
    result = 0
    for i in range(min(len(vector_left.components), len(vector_right.components))):
        result += (vector_left.components[i] * vector_right.components[i])
    return result
