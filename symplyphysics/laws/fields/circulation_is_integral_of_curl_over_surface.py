from typing import Sequence
from sympy import (Expr, Eq, Integral, simplify, Symbol as SymSymbol)
from symplyphysics import (Vector, dot_vectors, Quantity, print_expression)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.geometry.normals import parametrized_surface_normal
from symplyphysics.core.geometry.parameters import is_parametrized_surface
from symplyphysics.laws.fields import flux_is_integral_across_surface as flux_def

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
## C is circulation
## F is vector field
## Surface is surface boundary enclosed inside the curve
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

# These are not physical symbols - SymPy 'symbols' is good enough.

circulation = SymSymbol("circulation")
# This is integrand of the integral over the surface
curl_dot_surface_element = SymSymbol("curl_dot_surface_element")
# Inner integral parameter and limits
parameter1 = SymSymbol("parameter1")
parameter1_from = SymSymbol("parameter1_from")
parameter1_to = SymSymbol("parameter1_to")
# Outer integral parameter and limits
parameter2 = SymSymbol("parameter2")
parameter2_from = SymSymbol("parameter2_from")
parameter2_to = SymSymbol("parameter2_to")

law = Eq(
    circulation,
    Integral(curl_dot_surface_element, (parameter1, parameter1_from, parameter1_to),
    (parameter2, parameter2_from, parameter2_to)))

# Verify that circulation over surface equals to flux of curl of field.

flux_field_integral = flux_def.law.subs({
    flux_def.field_dot_surface_element: curl_dot_surface_element,
    flux_def.parameter1: parameter1,
    flux_def.parameter1_from: parameter1_from,
    flux_def.parameter1_to: parameter1_to,
    flux_def.parameter2: parameter2,
    flux_def.parameter2_from: parameter2_from,
    flux_def.parameter2_to: parameter2_to,
})
assert flux_field_integral.rhs == law.rhs


def print_law() -> str:
    return print_expression(law)


# Calculate SurfaceIntegral integrand, which is Dot(Curl(Field), dS)
def _calculate_curl_dot_surface_element(field_: VectorField, surface_: Sequence[Expr]) -> Expr:
    field_rotor_vector_field = curl_operator(field_)
    field_rotor_applied = field_rotor_vector_field.apply(surface_)
    surface_vector = Vector(surface_, field_.coordinate_system)
    surface_element_vector = parametrized_surface_normal(surface_vector, parameter1, parameter2)
    return dot_vectors(field_rotor_applied, surface_element_vector)


# field_ should be VectorField
# surface_ should be array with projections to coordinates, eg [parameter1 * cos(parameter2), parameter1 * sin(parameter2)]
def calculate_circulation(field_: VectorField, surface_: Sequence[Expr],
    parameter1_limits: tuple[ScalarValue, ScalarValue], parameter2_limits: tuple[ScalarValue,
    ScalarValue]) -> Quantity:
    if not is_parametrized_surface(surface_, parameter1, parameter2):
        raise ValueError("Trajectory should be parametrized by both parameter1 and parameter2")
    curl_dot_surface_element_value = _calculate_curl_dot_surface_element(field_, surface_)
    circulation_value = law.rhs.subs({
        curl_dot_surface_element: curl_dot_surface_element_value,
        parameter1_from: parameter1_limits[0],
        parameter1_to: parameter1_limits[1],
        parameter2_from: parameter2_limits[0],
        parameter2_to: parameter2_limits[1],
    }).doit()
    return Quantity(simplify(circulation_value))
