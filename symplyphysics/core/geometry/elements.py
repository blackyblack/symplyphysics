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


def volume_element_magnitude(volume: Vector) -> ScalarValue:
    if volume.coordinate_system == CoordinateSystem.System.CARTESIAN:
        x = volume.coordinate_system.coord_system.base_scalars()[0]
        y = volume.coordinate_system.coord_system.base_scalars()[1]
        z = volume.coordinate_system.coord_system.base_scalars()[2]
        return curve_element_magnitude(volume, x) * curve_element_magnitude(volume,
            y) * curve_element_magnitude(volume, z)
    if volume.coordinate_system == CoordinateSystem.System.CYLINDRICAL:
        r = volume.coordinate_system.coord_system.base_scalars()[0]
        theta = volume.coordinate_system.coord_system.base_scalars()[1]
        z = volume.coordinate_system.coord_system.base_scalars()[2]
        return r * curve_element_magnitude(volume, r) * curve_element_magnitude(volume,
            theta) * curve_element_magnitude(volume, z)
    if volume.coordinate_system == CoordinateSystem.System.SPHERICAL:
        r = volume.coordinate_system.coord_system.base_scalars()[0]
        theta = volume.coordinate_system.coord_system.base_scalars()[1]
        phi = volume.coordinate_system.coord_system.base_scalars()[2]
        return r**2 * sin(phi) * curve_element_magnitude(volume, r) * curve_element_magnitude(
            volume, theta) * curve_element_magnitude(volume, phi)
    raise ValueError(f"Unsupported coordinate system: {volume.coordinate_system}")
