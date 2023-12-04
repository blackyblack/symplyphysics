from typing import Sequence
from sympy import (Expr, Symbol as SymSymbol)
from symplyphysics import Quantity
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.analysis import flux_across_curve
from symplyphysics.core.fields.vector_field import VectorField

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
## dl is curve element magnitude
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

# This is how we parametrize the curve
parameter = SymSymbol("parameter")


def flux_law(field: VectorField, trajectory: Sequence[Expr], parameter_from: ScalarValue,
    parameter_to: ScalarValue) -> ScalarValue:
    return flux_across_curve(field, trajectory, (parameter, parameter_from, parameter_to))


# trajectory should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
# trajectory and field should be 2-dimensional, on XY plane
def calculate_flux(field: VectorField, trajectory: Sequence[Expr],
    parameter_limits: tuple[ScalarValue, ScalarValue]) -> Quantity:
    (parameter_from, parameter_to) = parameter_limits
    result_expr = flux_law(field, trajectory, parameter_from, parameter_to)
    return Quantity(result_expr)
