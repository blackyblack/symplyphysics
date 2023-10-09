from sympy import Expr, diff

from symplyphysics.core.geometry.elements import curve_element
from ..vectors.vectors import Vector
from ..vectors.arithmetics import cross_cartesian_vectors


# Curve normal is orthogonal to its tangent.
# Curve should be on XY plane.
def curve_normal(trajectory: Vector, parameter: Expr) -> Vector:
    curve_element_vector = curve_element(trajectory, parameter)
    k_vector = Vector([0, 0, 1], trajectory.coordinate_system)
    return cross_cartesian_vectors(curve_element_vector, k_vector)


# Calculate surface normal, which is Cross(Derivative(Surface, x), Derivative(Surface, y)) for Surface
# parametrized with 2 parameters
def parametrized_surface_normal(surface: Vector, parameter1: Expr, parameter2: Expr) -> Vector:
    surface_sympy_vector = surface.to_sympy_vector()
    surface_element_x = diff(surface_sympy_vector, parameter1)
    surface_element_vector_x = Vector.from_sympy_vector(surface_element_x,
        surface.coordinate_system)
    surface_element_y = diff(surface_sympy_vector, parameter2)
    surface_element_vector_y = Vector.from_sympy_vector(surface_element_y,
        surface.coordinate_system)
    return cross_cartesian_vectors(surface_element_vector_x, surface_element_vector_y)
