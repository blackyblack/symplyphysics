from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_vector_input, validate_vector_output,
    expr_to_vector_of_quantities
)
from symplyphysics.core.vectors.vectors import Vector, sympy_vector_from_vector

# Description
## Deformed sprign is about to return back to it's undeformed state and responds with some force. Law is:
## F = -kx, where
## F is force of response vector,
## k is elastic coefficient,
## x is vector of deformation.

# Condition
## Deformation is elactic (reversible).

response_force, elastic_coefficient, deformation = symbols('response_force elastic_coefficient deformation')
law = Eq(response_force, -1 * elastic_coefficient * deformation)

def print():
    return pretty(law, use_unicode=False)


@validate_input(coefficient_=units.force / units.length)
@validate_vector_input(deformation_=units.length)
@validate_vector_output(units.force)
def calculate_force(coefficient_: Quantity, deformation_: Vector) -> Vector:
    result_force_expr = solve(law, response_force, dict=True)[0][response_force]
    sympy_vector_deformation = sympy_vector_from_vector(deformation_)
    result_expr = result_force_expr.subs({elastic_coefficient: coefficient_, deformation: sympy_vector_deformation})
    return expr_to_vector_of_quantities(result_expr, "result_force", deformation_.coordinate_system)
    