from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.core.vectors.vectors import extended_express
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from sympy.vector import CoordSys3D

# Description
## Deformed sprign is about to return back to it's undeformated state and responds with some force. Law is:
## F = -kx, where
## F is force of response vector,
## k is elastic cpefficient,
## x is vector of deformation.

# Condition
## Deformation is elactic (reversible).

response_force, elastic_coefficient, deformation = symbols('response_force elastic_coefficient deformation')
law = Eq(response_force, scale_vector(-elastic_coefficient, deformation))

def print():
    return pretty(law, use_unicode=False)


@validate_input(coefficient_=units.force / units.length, deformation_=units.length, deformation_angle_=angle_type)
@validate_output(units.force)
def calculate_force(coefficient_: Quantity, deformation_: Quantity, deformation_angle_: Quantity) -> Quantity:
    result_force_expr = solve(law, response_force, dict=True)[0][response_force]
    Ca = CoordSys3D("Ca")
    Cy = Ca.create_new("Cy", transformation="cylindrical")
    #HACK: sympy angles are always in radians
    deformation_angle_radians = deformation_angle_.scale_factor
    
    # Convert to Cartesian coordinates as nothing works properly for Cylindrical coordinates
    def_vector = deformation_ * Cy.i + deformation_angle_radians * Cy.j    
    def_vector_cartesian = extended_express(def_vector, Ca)    
    result_expr = result_force_expr.subs({elastic_coefficient: coefficient_, deformation: def_vector_cartesian}).doit()
    return expr_to_quantity(result_expr, "force_vector")
    