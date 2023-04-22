from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    Vector,
    sympy_vector_from_vector,
    validate_input_symbols,
)
from symplyphysics.core.expr_to_quantity import expr_to_vector_of_quantities
from symplyphysics.core.quantity_decorator import validate_vector_input, validate_vector_output

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


def print() -> str:
    return print_expression(law)


@validate_input_symbols(coefficient_=elastic_coefficient)
@validate_vector_input(deformation_=deformation.dimension)
@validate_vector_output(response_force.dimension)
def calculate_force(coefficient_: Quantity, deformation_: Vector) -> Vector:
    result_force_expr = solve(law, response_force, dict=True)[0][response_force]
    sympy_vector_deformation = sympy_vector_from_vector(deformation_)
    result_expr = result_force_expr.subs({
        elastic_coefficient: coefficient_,
        deformation: sympy_vector_deformation
    })
    #TODO: think about some better solution to processing vectors
    return expr_to_vector_of_quantities(result_expr, deformation_.coordinate_system)
