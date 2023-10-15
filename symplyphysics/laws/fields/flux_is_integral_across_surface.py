from typing import Sequence
from sympy import (Expr, Symbol as SymSymbol)
from symplyphysics import Quantity
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.analysis import flux_across_surface
from symplyphysics.core.fields.vector_field import VectorField

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

# Inner integral parameter
parameter1 = SymSymbol("parameter1")
# Outer integral parameter
parameter2 = SymSymbol("parameter2")


def flux_law(field: VectorField, trajectory: Sequence[Expr],
    parameter1_limits: tuple[ScalarValue, ScalarValue],
    parameter2_limits: tuple[ScalarValue, ScalarValue]) -> ScalarValue:
    return flux_across_surface(field, trajectory, (parameter1, parameter1_limits[0], parameter1_limits[1]),
        (parameter2, parameter2_limits[0], parameter2_limits[1]))


# surface should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
def calculate_flux(field: VectorField, surface: Sequence[Expr],
    parameter1_limits: tuple[ScalarValue, ScalarValue], parameter2_limits: tuple[ScalarValue,
    ScalarValue]) -> Quantity:
    result_expr = flux_law(field, surface, parameter1_limits, parameter2_limits)
    return Quantity(result_expr)
