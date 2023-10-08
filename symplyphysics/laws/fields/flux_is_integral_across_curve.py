from typing import Sequence
from sympy import (Expr, Integral, Eq, simplify, Symbol as SymSymbol, sympify)
from symplyphysics import (cross_cartesian_vectors, print_expression, Quantity, Vector, dot_vectors,
    vector_unit, vector_magnitude)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.geometry.elements import curve_element

# Description
## Flux is defined as the amount of "stuff" going through a curve or a surface
## and we can get the flux at a particular point by taking the force and seeing
## how much of the force is perpendicular to the curve.
## Curve is parametrized counter-clockwise. Positive field flux assumes the direction of the field is outwards
## of the curve.

# Definition
## Flux = CurveIntegral(dot(F, n) * dl, Curve)
# Where:
## Flux is flux
## F is vector field
## n (norm) is outward unit normal vector of the curve
## dl is curve unit tangent vector magnitude
## dot is dot product

# Conditions
## - Curve is a function of a single parameter (eg y(x) = x**2), or parametrized with a single parameter
## (eg x(t) = cos(t), y(t) = sin(t))
## - Curve is smooth and continuous
## - Curve is positively oriented
## - Curve is defined for a plane XY (2-dimensional coordinate system)
## - Field is defined for a plane XY (2-dimensional coordinate system)

# NOTE:
## 3-dimensional curves, eg helix, are not supported because their norm is not well defined.

# These are not physical symbols - SymPy 'Symbol' is good enough.

flux = SymSymbol("flux")
curve_element_magnitude = SymSymbol("curve_element_magnitude")
parameter = SymSymbol("parameter")
parameter_from = SymSymbol("parameter_from")
parameter_to = SymSymbol("parameter_to")
field_dot_norm = SymSymbol("field_dot_norm")

law = Eq(
    flux,
    Integral(field_dot_norm * curve_element_magnitude, (parameter, parameter_from, parameter_to)))


def print_law() -> str:
    return print_expression(law)


def _unit_norm_vector(trajectory_: Vector):
    curve_element_vector = curve_element(trajectory_, parameter)
    unit_tangent_vector = vector_unit(curve_element_vector)
    # Use cross product with k vector to assert that norm vector is on
    # the XY plane
    k_vector = Vector([0, 0, 1], trajectory_.coordinate_system)
    return cross_cartesian_vectors(unit_tangent_vector, k_vector)


def _calculate_dot_norm(field_: VectorField, trajectory_: Vector) -> Expr:
    norm_vector = _unit_norm_vector(trajectory_)
    trajectory_components = [sympify(c) for c in trajectory_.components]
    field_applied = field_.apply(trajectory_components)
    return dot_vectors(field_applied, norm_vector)


def _curve_element_magnitude(trajectory_: Vector) -> Expr:
    curve_element_vector = curve_element(trajectory_, parameter)
    return vector_magnitude(curve_element_vector)


# trajectory_ should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
# trajectory_ and field_ should be 2-dimensional, on XY plane
def calculate_flux(field_: VectorField, trajectory_: Sequence[Expr],
    parameter_limits: tuple[ScalarValue, ScalarValue]) -> Quantity:
    if len(trajectory_) > 2:
        raise ValueError(f"Trajectory should have at most 2 components, got {len(trajectory_)}")
    trajectory_vector = Vector(trajectory_, field_.coordinate_system)
    field_dot_norm_value = _calculate_dot_norm(field_, trajectory_vector)
    curve_element_magnitude_value = _curve_element_magnitude(trajectory_vector)
    flux_value = law.rhs.subs({
        field_dot_norm: field_dot_norm_value,
        curve_element_magnitude: curve_element_magnitude_value,
        parameter_from: parameter_limits[0],
        parameter_to: parameter_limits[1],
    }).doit()
    return Quantity(simplify(flux_value))
