from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, Vector,
    sympy_vector_from_vector, validate_input, validate_output)
from symplyphysics.core.vectors.vectors import expr_to_vector

# Description
## Deformed sprign is about to return back to it's undeformed state and responds with some force. Law is:
## F = -kx, where
## F is force of response vector,
## k is elastic coefficient,
## x is vector of deformation.

# Condition
## Deformation is elactic (reversible).

response_force = Symbol("response_force", units.force)
elastic_coefficient = Symbol("elastic_coefficient", units.force / units.length)
deformation = Symbol("deformation", units.length)

law = Eq(response_force, -1 * elastic_coefficient * deformation)


def print_law() -> str:
    return print_expression(law)


@validate_input(coefficient_=elastic_coefficient, deformation_=deformation)
@validate_output(response_force)
def calculate_force(coefficient_: Quantity, deformation_: Vector) -> Vector:
    result_force_expr = solve(law, response_force, dict=True)[0][response_force]
    sympy_vector_deformation = sympy_vector_from_vector(deformation_)
    result_expr = result_force_expr.subs({
        elastic_coefficient: coefficient_,
        deformation: sympy_vector_deformation
    })
    return expr_to_vector(result_expr, deformation_.coordinate_system)
