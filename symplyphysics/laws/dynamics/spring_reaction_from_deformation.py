from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units,
    validate_vector_input, expr_to_vector_of_quantities, validate_vector_output
)
from symplyphysics.core.quantity_decorator import validate_input_symbols
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable
from symplyphysics.core.vectors.vectors import Vector, sympy_vector_from_vector

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

def print(expr: Expr) -> str:
    symbols = [response_force, elastic_coefficient, deformation]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(coefficient_=elastic_coefficient)
@validate_vector_input(deformation_=deformation.dimension)
@validate_vector_output(response_force.dimension)
def calculate_force(coefficient_: Quantity, deformation_: Vector) -> Vector:
    result_force_expr = solve(law, response_force, dict=True)[0][response_force]
    sympy_vector_deformation = sympy_vector_from_vector(deformation_)
    result_expr = result_force_expr.subs({elastic_coefficient: coefficient_, deformation: sympy_vector_deformation})
    #TODO: think about some better solution to processing vectors
    return expr_to_vector_of_quantities(result_expr, "response_force", deformation_.coordinate_system)
