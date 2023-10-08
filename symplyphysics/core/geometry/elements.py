from sympy import Expr, diff
from ..vectors.vectors import Vector
from ..vectors.arithmetics import cross_cartesian_vectors


# Curve element is its tangent vector
def curve_element(trajectory: Vector, parameter: Expr) -> Vector:
    trajectory_sympy_vector = trajectory.to_sympy_vector()
    trajectory_element_sympy_vector = diff(trajectory_sympy_vector, parameter)
    return Vector.from_sympy_vector(trajectory_element_sympy_vector, trajectory.coordinate_system)


# Calculate vector surface element, which is Cross(Derivative(Surface, x), Derivative(Surface, y)) for Surface
# parametrized with 2 parameters
def parametrized_surface_element(surface: Vector, parameter1: Expr, parameter2: Expr) -> Vector:
    surface_sympy_vector = surface.to_sympy_vector()
    surface_element_x = diff(surface_sympy_vector, parameter1)
    surface_element_vector_x = Vector.from_sympy_vector(surface_element_x,
        surface.coordinate_system)
    surface_element_y = diff(surface_sympy_vector, parameter2)
    surface_element_vector_y = Vector.from_sympy_vector(surface_element_y,
        surface.coordinate_system)
    return cross_cartesian_vectors(surface_element_vector_x, surface_element_vector_y)
