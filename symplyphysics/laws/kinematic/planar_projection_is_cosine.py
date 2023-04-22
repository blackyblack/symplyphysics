import numbers
from sympy import (Eq, solve, symbols, cos)
from symplyphysics import (
    expr_to_quantity, Quantity, Symbol, print_expression, angle_type,
    validate_input_symbols,
)
from symplyphysics.core.quantity_decorator import validate_output_same

# Description
## Most of cases might be represented in 2-dimensional space with two orthogonal axis - vertical Y and horizontal X. Any vector in this space (velocity, force etc) can be easily
## transformed to sum of two orthogonal vectors - projections of this vector to axis with help of angle between the vector and the axis.
## Law: Vaxis = V * cos(alpha)
## Where:
## Vaxis is the projection of vector V to the axis
## V is length of projected vector
## alpha is the angle between vector and axis

vector_angle = Symbol("vector_angle", angle_type)
vector_length = symbols("vector_length")
projection = symbols("projection")

law = Eq(projection, vector_length * cos(vector_angle))

def print() -> str:
    return print_expression(law)

@validate_input_symbols(angle_=vector_angle)
@validate_output_same("vector_length_")
def calculate_projection(vector_length_: Quantity, angle_: Quantity | float) -> Quantity:
    result_projection_expr = solve(law, projection, dict=True)[0][projection]
    #HACK: sympy angles are always in radians
    angle_radians = angle_ if isinstance(angle_, numbers.Number) else angle_.scale_factor
    result_expr = result_projection_expr.subs({vector_length: vector_length_, vector_angle: angle_radians})
    return expr_to_quantity(result_expr)
