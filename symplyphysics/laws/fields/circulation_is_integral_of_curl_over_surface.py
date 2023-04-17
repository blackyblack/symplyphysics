from typing import List
from sympy import Integral, Derivative
from sympy.vector import Dot, Curl, Cross
from symplyphysics import (
    symbols, Eq, pretty, VectorField, simplify
)
from symplyphysics.core.fields.vector_field import field_from_sympy_vector, sympy_vector_from_field
from symplyphysics.core.vectors.vectors import Vector, sympy_vector_from_vector

# Description
## Circulation of the field along the closed curve is flow of the rotor (or curl) of this field 
## through any area surrounded by this curve.
## CurveIntegral(F * dl, Curve) = SurfaceIntegral(Curl(F) * dS, Surface), where
## S is area surrounded by Curve.
## Potential field is the field with zero rotor. Also potential field is called irrotational field.
## Work to move the object along the closed curve in the potential field is zero.

# Definition
## C = SurfaceIntegral(Curl(F) * dS, Surface), where
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

circulation, field, surface = symbols("circulation field surface")
# field rotor should be evaluated and applied before using in circulation definition
field_rotor = symbols("field_rotor")
# surface_element (dS) is surface derivative by two parameters
surface_element, surface_element_by_parameter1, surface_element_by_parameter2 = symbols("surface_element surface_element_by_parameter1 surface_element_by_parameter2")
parameter1, parameter1_from, parameter1_to = symbols("parameter1 parameter1_from parameter1_to")
parameter2, parameter2_from, parameter2_to = symbols("parameter2 parameter2_from parameter2_to")

# field_rotor, surface_element_by_parameter1, surface_element_by_parameter2, surface_element - should be evaluated before passing to definition
# see calculate_circulation() for an example
field_rotor_definition = Eq(field_rotor, Curl(field))
surface_element_by_parameter1_definition = Eq(surface_element_by_parameter1, Derivative(surface, parameter1))
surface_element_by_parameter2_definition = Eq(surface_element_by_parameter2, Derivative(surface, parameter2))
surface_element_definition = Eq(surface_element, Cross(surface_element_by_parameter1, surface_element_by_parameter2), evaluate=False)
definition = Eq(circulation, Integral(Dot(field_rotor, surface_element), (parameter1, parameter1_from, parameter1_to), (parameter2, parameter2_from, parameter2_to)))

def print():
    return pretty(definition, use_unicode=False)

# field_ should be VectorField
# surface_ should be array with projections to coordinates, eg [parameter1 * cos(parameter2), parameter1 * sin(parameter2)]
def calculate_circulation(
    field_: VectorField,
    surface_: List,
    parameter1_from_,
    parameter1_to_,
    parameter2_from_,
    parameter2_to_):

    field_space = sympy_vector_from_field(field_)
    field_rotor_sympy = field_rotor_definition.rhs.subs(field, field_space).doit()
    field_rotor_applied = field_from_sympy_vector(field_rotor_sympy, field_.coordinate_system)
    field_applied = field_rotor_applied.apply(surface_)
    field_as_vector = sympy_vector_from_vector(field_applied)
    surface_sympy_vector = sympy_vector_from_vector(Vector(surface_, field_.coordinate_system))
    surface_element_x = surface_element_by_parameter1_definition.rhs.subs(surface, surface_sympy_vector).doit()
    surface_element_y = surface_element_by_parameter2_definition.rhs.subs(surface, surface_sympy_vector).doit()
    surface_element_result = surface_element_definition.rhs.subs({surface_element_by_parameter1: surface_element_x, surface_element_by_parameter2: surface_element_y}).doit()
    result_expr = definition.rhs.subs({field_rotor: field_as_vector, surface_element: surface_element_result, parameter1_from: parameter1_from_, parameter1_to: parameter1_to_, parameter2_from: parameter2_from_, parameter2_to: parameter2_to_}).doit()
    return simplify(result_expr)
