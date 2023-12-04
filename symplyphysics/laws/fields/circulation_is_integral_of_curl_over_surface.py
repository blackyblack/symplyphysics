from typing import Sequence
from sympy import (Expr, Symbol as SymSymbol)
from symplyphysics import Quantity
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.analysis import circulation_along_surface_boundary
from symplyphysics.core.fields.vector_field import VectorField

# Description
## Circulation of the field along the closed curve is flow of the rotor (or curl) of this field
## through any area surrounded by this curve.
## See [circulation definition](..\..\definitions\circulation_is_integral_along_curve.py) for circulation
## as curvilinear integral formula.
## Surface integral is derived from linear integral by using [Stoke's theorem](https://en.wikipedia.org/wiki/Stokes%27_theorem).
## Potential field is the field with zero rotor. Also potential field is called irrotational field.
## Work to move the object along the closed curve in the potential field is zero.
## In other words, the line integral of a vector field over a loop is equal to the flux of its curl through the enclosed surface.
## See [flux definition](./flux_is_integral_of_across_surface.py) for more information on vector field flux.

# Law:
## C = SurfaceIntegral(dot(Curl(F), dS), Surface)
## Where:
## C is circulation along Curve
## F is vector field
## Surface is surface enclosed inside the Curve
## dS is vector surface element, equals to n * dA, where:
## - n is unit normal vector of the surface
## - dA is infinitesimal element of surface area
## dot is dot product
## Curl is rotor (or curl) operator

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


def circulation_law(field: VectorField, trajectory: Sequence[Expr],
    parameter1_limits: tuple[ScalarValue, ScalarValue], parameter2_limits: tuple[ScalarValue,
    ScalarValue]) -> ScalarValue:
    return circulation_along_surface_boundary(field, trajectory,
        (parameter1, parameter1_limits[0], parameter1_limits[1]),
        (parameter2, parameter2_limits[0], parameter2_limits[1]))


# surface should be array with projections to coordinates, eg [parameter1 * cos(parameter2), parameter1 * sin(parameter2)]
def calculate_circulation(field: VectorField, surface: Sequence[Expr],
    parameter1_limits: tuple[ScalarValue, ScalarValue], parameter2_limits: tuple[ScalarValue,
    ScalarValue]) -> Quantity:
    result_expr = circulation_law(field, surface, parameter1_limits, parameter2_limits)
    return Quantity(result_expr)
