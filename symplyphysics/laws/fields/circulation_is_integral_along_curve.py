from typing import Sequence
from sympy import (Expr, Symbol as SymSymbol)
from symplyphysics import Quantity
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.analysis import circulation_along_curve
from symplyphysics.core.fields.vector_field import VectorField

# Description
## Field circulation along closed curve is curvilinear integral of field by this curve.
## A space curve is parametrized by a vector-valued function that depends upon a single parameter t that varies over some interval.
## We shall always assume that x(t) is continuously differentiable. The curve is smooth provided its tangent
## vector is continuous and everywhere nonzero.
## Physically, we can think of a curve as the trajectory described by a particle moving in
## space. At each time t, the tangent vector x'(t) represents the instantaneous velocity of the
## particle. Thus, as long as the particle moves with nonzero speed, its trajectory is necessarily a smooth curve.
## When curve is not smooth (eg trajectory is a square) one should represent this curve as a sum of its smooth parts
## (eg sum of sides of a square).
## Curve is parametrized counter-clockwise. Positive field circulation assumes the direction of the field is the same as
## the curve orientation.

# Definition
## C = CurveIntegral(dot(F, dr), Curve)
# Where:
## C is circulation
## F is vector field
## dr is curve element (tangent) vector
## dot is dot product

# Conditions
## - Curve is a function of a single parameter (eg y(x) = x**2), or parametrized with a single parameter
## (eg x(t) = cos(t), y(t) = sin(t))
## - Curve is smooth and continuous
## - Curve is positively oriented

# This is how we parametrize the curve
parameter = SymSymbol("parameter")


def circulation_law(field: VectorField, trajectory: Sequence[Expr], parameter_from: ScalarValue,
    parameter_to: ScalarValue) -> ScalarValue:
    return circulation_along_curve(field, trajectory, (parameter, parameter_from, parameter_to))


# trajectory should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
def calculate_circulation(field: VectorField, trajectory: Sequence[Expr],
    parameter_limits: tuple[ScalarValue, ScalarValue]) -> Quantity:
    (parameter_from, parameter_to) = parameter_limits
    result_expr = circulation_law(field, trajectory, parameter_from, parameter_to)
    return Quantity(result_expr)
