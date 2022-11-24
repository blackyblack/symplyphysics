from sympy import Integral, Derivative
from sympy.vector import Dot, Curl, Cross
from symplyphysics import (
    symbols, Eq, pretty, simplify, apply_field
)

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

circulation, field, surface = symbols('circulation field surface')
# field rotor should be evaluated and applied before using in circulation definition
field_rotor = symbols('field_rotor')
# surface_element (dS) is surface derivative by two parameters
surface_element, surface_element_by_parameter1, surface_element_by_parameter2 = symbols('surface_element surface_element_by_parameter1 surface_element_by_parameter2')
parameter1, parameter1_from, parameter1_to = symbols('parameter1 parameter1_from parameter1_to')
parameter2, parameter2_from, parameter2_to = symbols('parameter2 parameter2_from parameter2_to')

# field_rotor, surface_element_by_parameter1, surface_element_by_parameter2, surface_element - should be evaluated before passing to definition
# see calculate_circulation() for an example
field_rotor_definition = Eq(field_rotor, Curl(field))
surface_element_by_parameter1_definition = Eq(surface_element_by_parameter1, Derivative(surface, parameter1))
surface_element_by_parameter2_definition = Eq(surface_element_by_parameter2, Derivative(surface, parameter2))
surface_element_definition = Eq(surface_element, Cross(surface_element_by_parameter1, surface_element_by_parameter2), evaluate=False)
definition = Eq(circulation, Integral(Dot(field_rotor, surface_element), (parameter1, parameter1_from, parameter1_to), (parameter2, parameter2_from, parameter2_to)))

def print():
    return pretty(definition, use_unicode=False)

# field_ should be instance of CoordSys3D, eg C.y * C.i + -1 * C.x * C.j
# surface_ should be instance of CoordSys3D, eg parameter1 * cos(parameter2) * C.i + parameter1 * sin(parameter2) * C.j
def calculate_circulation(
    field_,
    surface_,
    parameter1_from_,
    parameter1_to_,
    parameter2_from_,
    parameter2_to_):

    field_rotor_applied = field_rotor_definition.rhs.subs(field, field_).doit()
    field_applied = apply_field(field_rotor_applied, surface_)
    surface_element_x = surface_element_by_parameter1_definition.rhs.subs(surface, surface_).doit()
    surface_element_y = surface_element_by_parameter2_definition.rhs.subs(surface, surface_).doit()
    surface_element_result = surface_element_definition.rhs.subs({surface_element_by_parameter1: surface_element_x, surface_element_by_parameter2: surface_element_y}).doit()
    result_expr = definition.rhs.subs({field_rotor: field_applied, surface_element: surface_element_result, parameter1_from: parameter1_from_, parameter1_to: parameter1_to_, parameter2_from: parameter2_from_, parameter2_to: parameter2_to_}).doit()
    # some expressions are invalid without simplifying them first
    return simplify(result_expr)
