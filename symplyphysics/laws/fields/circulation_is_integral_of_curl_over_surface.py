from typing import Sequence
from sympy import (Expr, Eq, Integral, Derivative, simplify, Symbol as SymSymbol)
from symplyphysics import (CoordinateSystem, Vector, cross_cartesian_vectors, dot_vectors, Quantity)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField

# Description
## Circulation of the field along the closed curve is flow of the rotor (or curl) of this field
## through any area surrounded by this curve.
## CurveIntegral(F * dl, Curve) = SurfaceIntegral(Curl(F) * dS, Surface), where
## S is area surrounded by Curve.
## Potential field is the field with zero rotor. Also potential field is called irrotational field.
## Work to move the object along the closed curve in the potential field is zero.

# Law:
## C = SurfaceIntegral(Curl(F) * dS, Surface)
## Where:
## C is circulation
## F is vector field
## S is surface boundary, equals to curve area
## dS is surface double derivative
## * is dot product
## Curl is rotor (or curl) operator

# Conditions
## - Field is smooth vector field in 3d space
## - Surface is smooth oriented surface in 3d space
## - Curve is smooth, continuous and closed

# These are not physical symbols - SymPy 'symbols' is good enough.

surface = SymSymbol("surface")

# surface_element (dS) is surface derivative by two parameters
# surface_partial_element is surface derivative by 1 coordinate
surface_partial_element = SymSymbol("surface_partial_element")
surface_element_parameter = SymSymbol("surface_element_parameter")

surface_partial_element_definition = Eq(surface_partial_element,
    Derivative(surface, surface_element_parameter))


# Calculate surface element, which is Cross(Derivative(Surface, x), Derivative(Surface, y))
def _surface_element(surface_: Sequence[Expr], parameters: Sequence[SymSymbol],
    coordinate_system: CoordinateSystem) -> Vector:
    assert len(parameters) == 2, "Surface element should have 2 parameters"
    surface_sympy_vector = Vector(surface_, coordinate_system).to_sympy_vector()
    surface_element_x = surface_partial_element_definition.rhs.subs(surface_element_parameter,
        parameters[0])
    surface_element_x = surface_element_x.subs(surface, surface_sympy_vector).doit()
    surface_element_y = surface_partial_element_definition.rhs.subs(surface_element_parameter,
        parameters[1])
    surface_element_y = surface_element_y.subs(surface, surface_sympy_vector).doit()
    surface_element_vector_x = Vector.from_sympy_vector(surface_element_x, coordinate_system)
    surface_element_vector_y = Vector.from_sympy_vector(surface_element_y, coordinate_system)
    return cross_cartesian_vectors(surface_element_vector_x, surface_element_vector_y)


# Calculate SurfaceIntegral integrand, which is Dot(Curl(Field), dS)
def integral_expression(field_: VectorField, surface_: Sequence[Expr],
    parameters: Sequence[SymSymbol]) -> Expr:
    surface_element_vector = _surface_element(surface_, parameters, field_.coordinate_system)
    field_rotor_vector_field = curl_operator(field_)
    field_rotor_applied = field_rotor_vector_field.apply(surface_)
    return dot_vectors(field_rotor_applied, surface_element_vector)


def circulation_law(integral_expression_: Expr, parameters: Sequence[SymSymbol],
    limits: Sequence[tuple[ScalarValue, ScalarValue]]) -> Expr:
    integral_limits = [(p, limit1, limit2) for (p, (limit1, limit2)) in zip(parameters, limits)]
    return Integral(integral_expression_, *integral_limits).doit()


# field_ should be VectorField
# surface_ should be array with projections to coordinates, eg [parameter1 * cos(parameter2), parameter1 * sin(parameter2)]
def calculate_circulation(field_: VectorField, surface_: Sequence[Expr],
    parameters: Sequence[SymSymbol], limits: Sequence[tuple[ScalarValue, ScalarValue]]) -> Quantity:
    integral_expression_value = integral_expression(field_, surface_, parameters)
    result_expr = circulation_law(integral_expression_value, parameters, limits)
    return Quantity(simplify(result_expr))
