from typing import Sequence
from sympy import (Basic, Expr, Eq, Integral, Derivative, simplify, Symbol as SymSymbol)
from symplyphysics import (CoordinateSystem, Vector, cross_cartesian_vectors, dot_vectors, Quantity,
    print_expression)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.operators import curl_operator
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
## - Surface is a function of at most two parameters (eg z(x, y) = sqrt(x**2 + y**2)), or parametrized with at most
##   two parameters (eg x(t1, t2) = t1 * cos(t2), y(t1, t2) = t1 * sin(t2), z(t1, t2) = t1)

# These are not physical symbols - SymPy 'symbols' is good enough.

surface = SymSymbol("surface")
circulation = SymSymbol("circulation")

# surface_element (dS) is surface derivative by two parameters
# surface_partial_element is surface derivative by 1 coordinate
surface_partial_element = SymSymbol("surface_partial_element")
surface_element_parameter = SymSymbol("surface_element_parameter")

# This is integrand of the integral over the surface
curl_dot_surface_element = SymSymbol("curl_dot_surface_element")
# This is inner integral of the double integral
curl_dot_surface_integral = SymSymbol("curl_dot_surface_integral")

# Inner integral parameter and limits
parameter1 = SymSymbol("parameter1")
parameter1_from = SymSymbol("parameter1_from")
parameter1_to = SymSymbol("parameter1_to")
# Outer integral parameter and limits
parameter2 = SymSymbol("parameter2")
parameter2_from = SymSymbol("parameter2_from")
parameter2_to = SymSymbol("parameter2_to")

# Surface partial derivative is derivative of the 'surface' over some parameter
surface_partial_derivative_definition = Eq(surface_partial_element,
    Derivative(surface, surface_element_parameter))

curl_dot_surface_integral_definition = Eq(
    curl_dot_surface_integral,
    Integral(curl_dot_surface_element, (parameter1, parameter1_from, parameter1_to)))

law = Eq(circulation,
    Integral(curl_dot_surface_integral, (parameter2, parameter2_from, parameter2_to)))


def print_law() -> str:
    return print_expression(law)


# Calculate vector surface element, which is Cross(Derivative(Surface, x), Derivative(Surface, y)) for Surface
# parametrized with 2 parameters
def _surface_element(surface_: Sequence[Expr], coordinate_system: CoordinateSystem) -> Vector:
    surface_vector = Vector(surface_, coordinate_system)
    surface_sympy_vector = surface_vector.to_sympy_vector()
    surface_element_x = surface_partial_derivative_definition.rhs.subs(
        surface_element_parameter, parameter1)
    surface_element_x = surface_element_x.subs(surface, surface_sympy_vector).doit()
    surface_element_vector_x = Vector.from_sympy_vector(surface_element_x, coordinate_system)
    surface_element_y = surface_partial_derivative_definition.rhs.subs(
        surface_element_parameter, parameter2)
    surface_element_y = surface_element_y.subs(surface, surface_sympy_vector).doit()
    surface_element_vector_y = Vector.from_sympy_vector(surface_element_y, coordinate_system)
    return cross_cartesian_vectors(surface_element_vector_x, surface_element_vector_y)


# Calculate SurfaceIntegral integrand, which is Dot(Curl(Field), dS)
def _calculate_curl_dot_surface_element(field_: VectorField, surface_: Sequence[Expr]) -> Expr:
    surface_element_vector = _surface_element(surface_, field_.coordinate_system)
    field_rotor_vector_field = curl_operator(field_)
    field_rotor_applied = field_rotor_vector_field.apply(surface_)
    return dot_vectors(field_rotor_applied, surface_element_vector)


# field_ should be VectorField
# surface_ should be array with projections to coordinates, eg [parameter1 * cos(parameter2), parameter1 * sin(parameter2)]
def calculate_circulation(field_: VectorField, surface_: Sequence[Expr],
    parameter1_limits: tuple[ScalarValue, ScalarValue],
    parameter2_limits: tuple[ScalarValue, ScalarValue]) -> Quantity:
    free_symbols: set[Basic] = set()
    for ss in surface_:
        free_symbols = free_symbols.union(ss.free_symbols)
    if (parameter1 not in free_symbols) or (parameter2 not in free_symbols):
        raise ValueError("Trajectory should be parametrized by both parameter1 and parameter2")

    curl_dot_surface_element_value = _calculate_curl_dot_surface_element(
        field_, surface_)
    curl_dot_surface_integral_value = curl_dot_surface_integral_definition.rhs.subs({
        curl_dot_surface_element: curl_dot_surface_element_value,
        parameter1_from: parameter1_limits[0],
        parameter1_to: parameter1_limits[1],
    }).doit()
    circulation_value = law.rhs.subs({
        curl_dot_surface_integral: curl_dot_surface_integral_value,
        parameter2_from: parameter2_limits[0],
        parameter2_to: parameter2_limits[1],
    }).doit()
    return Quantity(simplify(circulation_value))
