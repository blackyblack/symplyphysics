from sympy import Expr, diff
from ..vectors.vectors import Vector


# Curve element is its tangent vector
def curve_element(trajectory: Vector, parameter: Expr) -> Vector:
    trajectory_sympy_vector = trajectory.to_sympy_vector()
    trajectory_element_sympy_vector = diff(trajectory_sympy_vector, parameter)
    return Vector.from_sympy_vector(trajectory_element_sympy_vector, trajectory.coordinate_system)
