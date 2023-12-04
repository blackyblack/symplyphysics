from sympy import Expr, diff, sin
from symplyphysics import CoordinateSystem, vector_magnitude
from ..dimensions import ScalarValue
from ..vectors.vectors import Vector


# Curve element is its tangent vector
def curve_element(trajectory: Vector, parameter: Expr) -> Vector:
    trajectory_sympy_vector = trajectory.to_sympy_vector()
    trajectory_element_sympy_vector = diff(trajectory_sympy_vector, parameter)
    return Vector.from_sympy_vector(trajectory_element_sympy_vector, trajectory.coordinate_system)


def curve_element_magnitude(trajectory: Vector, parameter: Expr) -> ScalarValue:
    return vector_magnitude(curve_element(trajectory, parameter))


# For non parametrized volumes
def volume_element_magnitude(coordinate_system: CoordinateSystem) -> ScalarValue:
    if coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        return 1
    if coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        r = coordinate_system.coord_system.base_scalars()[0]
        return r
    if coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        r = coordinate_system.coord_system.base_scalars()[0]
        phi = coordinate_system.coord_system.base_scalars()[2]
        return r**2 * sin(phi)
    raise ValueError(f"Unsupported coordinate system: {coordinate_system}")
