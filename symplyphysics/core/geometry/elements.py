from sympy import Expr, diff, sin, S
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.vectors.arithmetics import vector_magnitude
from ..vectors.vectors import Vector


# Curve element is its tangent vector
def parametrized_curve_element(trajectory: Vector, parameter: Expr) -> Vector:
    trajectory_sympy_vector = trajectory.to_sympy_vector()
    trajectory_element_sympy_vector = diff(trajectory_sympy_vector, parameter)
    return Vector.from_sympy_vector(trajectory_element_sympy_vector, trajectory.coordinate_system)


def parametrized_curve_element_magnitude(trajectory: Vector, parameter: Expr) -> Expr:
    return vector_magnitude(parametrized_curve_element(trajectory, parameter))


# For non parametrized volumes
def volume_element_magnitude(coordinate_system: CoordinateSystem) -> Expr:
    if coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        return S.One
    if coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        r = coordinate_system.coord_system.base_scalars()[0]
        return r
    if coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        r = coordinate_system.coord_system.base_scalars()[0]
        phi = coordinate_system.coord_system.base_scalars()[2]
        return r**2 * sin(phi)
    raise ValueError(f"Unsupported coordinate system: {coordinate_system}")
