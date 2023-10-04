from typing import Sequence
from sympy import (Expr, Integral, Derivative, Eq, simplify, Symbol as SymSymbol, sympify)
from symplyphysics import (cross_cartesian_vectors, print_expression, Quantity, Vector, dot_vectors,
    vector_unit, vector_magnitude)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.vector_field import VectorField

# Description
## Flux is defined as the amount of "stuff" going through a curve or a surface
## and we can get the flux at a particular point by taking the force and seeing
## how much of the force is perpendicular to the surface.
## Surface is parametrized counter-clockwise. Positive field flux assumes the direction of the field is outwards
## of the surface.

# Definition
## Flux = SurfaceIntegral(dot(F, n) * dS, Surface)
# Where:
## Flux is flux
## F is vector field
## n (norm) is outward unit normal vector of the surface
## dl is curve unit tangent vector magnitude
## dot is dot product

# Conditions
## - Curve is a function of a single parameter (eg y(x) = x**2), or parametrized with a single parameter
## (eg x(t) = cos(t), y(t) = sin(t))
## - Curve is smooth and continuous
## - Curve is positively oriented
## - Field is defined for a plane (2-dimensional coordinate system)

# These are not physical symbols - SymPy 'Symbol' is good enough.

flux = SymSymbol("flux")
# trajectory is a function of the moving particle
trajectory = SymSymbol("trajectory")
# trajectory_element (dl) is trajectory derivative by parameter
# parameter is an argument of trajectory function, eg x coordinate
trajectory_element = SymSymbol("trajectory_element")
trajectory_element_magnitude = SymSymbol("trajectory_element_magnitude")
parameter = SymSymbol("parameter")
parameter_from = SymSymbol("parameter_from")
parameter_to = SymSymbol("parameter_to")
field_dot_norm = SymSymbol("field_dot_norm")

trajectory_element_definition = Eq(trajectory_element, Derivative(trajectory, parameter))
law = Eq(
    flux,
    Integral(field_dot_norm * trajectory_element_magnitude,
    (parameter, parameter_from, parameter_to)))


def print_law() -> str:
    return print_expression(law)


def _unit_norm_vector(trajectory_: Vector):
    trajectory_sympy_vector = trajectory_.to_sympy_vector()
    tangent_sympy_vector = trajectory_element_definition.rhs.subs(trajectory,
        trajectory_sympy_vector).doit()
    tangent_vector = Vector.from_sympy_vector(tangent_sympy_vector, trajectory_.coordinate_system)
    unit_tangent_vector = vector_unit(tangent_vector)
    # Use cross product with k vector to assert that norm vector is not directed
    # upwards (in the same direction as k)
    k_vector = Vector([0, 0, 1], trajectory_.coordinate_system)
    return cross_cartesian_vectors(unit_tangent_vector, k_vector)


def _calculate_dot_norm(field_: VectorField, trajectory_: Vector) -> Expr:
    norm_vector = _unit_norm_vector(trajectory_)
    trajectory_components = [sympify(c) for c in trajectory_.components]
    field_applied = field_.apply(trajectory_components)
    return dot_vectors(field_applied, norm_vector)


def _trajectory_element_magnitude(trajectory_: Vector) -> Expr:
    trajectory_sympy_vector = trajectory_.to_sympy_vector()
    trajectory_element_sympy_vector = trajectory_element_definition.rhs.subs(
        trajectory, trajectory_sympy_vector).doit()
    trajectory_element_vector = Vector.from_sympy_vector(trajectory_element_sympy_vector,
        trajectory_.coordinate_system)
    return vector_magnitude(trajectory_element_vector)


# trajectory_ should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
def calculate_flux(field_: VectorField, trajectory_: Sequence[Expr],
    parameter_limits: tuple[ScalarValue, ScalarValue]) -> Quantity:
    # trajectory_ and field_ should be 2-dimensional
    if len(trajectory_) > 2:
        raise ValueError(f"Trajectory with 2 and more components is not supported,"
            f" got {len(trajectory_)}")
    trajectory_vector = Vector(trajectory_, field_.coordinate_system)
    field_dot_norm_value = _calculate_dot_norm(field_, trajectory_vector)
    trajectory_element_magnitude_value = _trajectory_element_magnitude(trajectory_vector)
    flux_value = law.rhs.subs({
        field_dot_norm: field_dot_norm_value,
        trajectory_element_magnitude: trajectory_element_magnitude_value,
        parameter_from: parameter_limits[0],
        parameter_to: parameter_limits[1],
    }).doit()
    return Quantity(simplify(flux_value))
