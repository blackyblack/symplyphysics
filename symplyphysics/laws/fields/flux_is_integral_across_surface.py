from typing import Sequence
from sympy import (Expr, Integral, Eq, simplify, Symbol as SymSymbol)
from symplyphysics import (print_expression, Quantity, Vector, dot_vectors)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.geometry.normals import parametrized_surface_normal
from symplyphysics.core.geometry.parameters import is_parametrized_surface

# Description
## Flux is defined as the amount of "stuff" going through a curve or a surface
## and we can get the flux at a particular point by taking the force and seeing
## how much of the force is perpendicular to the surface.
## Surface is parametrized counter-clockwise. Positive field flux assumes the direction of the field is outwards
## of the surface.
## See [flux across curve definition](.\flux_is_integral_across_curve.py) for flux as curvilinear
## integral formula.

# Definition
## Flux = SurfaceIntegral(dot(F, dS), Surface)
# Where:
## Flux is flux
## F is vector field
## dS is vector surface element, equals to n * dA, where:
## - n is unit normal vector of the surface
## - dA is infinitesimal element of surface area
## dot is dot product

# Conditions
## - Field is smooth vector field in 3d space
## - Surface is smooth positively oriented surface in 3d space
## - Surface enclosing curve is smooth, continuous and closed
## - Surface is a function of two parameters (eg z(x, y) = sqrt(x**2 + y**2)), or parametrized with
##   two parameters (eg x(t1, t2) = t1 * cos(t2), y(t1, t2) = t1 * sin(t2), z(t1, t2) = t1)

# These are not physical symbols - SymPy 'Symbol' is good enough.

flux = SymSymbol("flux")
# This is integrand of the integral over the surface
field_dot_surface_element = SymSymbol("field_dot_surface_element")
# Inner integral parameter and limits
parameter1 = SymSymbol("parameter1")
parameter1_from = SymSymbol("parameter1_from")
parameter1_to = SymSymbol("parameter1_to")
# Outer integral parameter and limits
parameter2 = SymSymbol("parameter2")
parameter2_from = SymSymbol("parameter2_from")
parameter2_to = SymSymbol("parameter2_to")

law = Eq(
    flux,
    Integral(field_dot_surface_element, (parameter1, parameter1_from, parameter1_to),
    (parameter2, parameter2_from, parameter2_to)))


def print_law() -> str:
    return print_expression(law)


# Calculate SurfaceIntegral integrand, which is Dot(Field, dS)
def _calculate_field_dot_surface_element(field_: VectorField, surface_: Sequence[Expr]) -> Expr:
    field_applied = field_.apply(surface_)
    surface_vector = Vector(surface_, field_.coordinate_system)
    surface_element_vector = parametrized_surface_normal(surface_vector, parameter1, parameter2)
    return dot_vectors(field_applied, surface_element_vector)


# trajectory_ should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
def calculate_flux(field_: VectorField, surface_: Sequence[Expr],
    parameter1_limits: tuple[ScalarValue, ScalarValue], parameter2_limits: tuple[ScalarValue,
    ScalarValue]) -> Quantity:
    if not is_parametrized_surface(surface_, parameter1, parameter2):
        raise ValueError("Trajectory should be parametrized by both parameter1 and parameter2")
    field_dot_surface_element_value = _calculate_field_dot_surface_element(field_, surface_)
    flux_value = law.rhs.subs({
        field_dot_surface_element: field_dot_surface_element_value,
        parameter1_from: parameter1_limits[0],
        parameter1_to: parameter1_limits[1],
        parameter2_from: parameter2_limits[0],
        parameter2_to: parameter2_limits[1],
    }).doit()
    return Quantity(simplify(flux_value))
