from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, cos, pi
)
#from sympy.physics.units.definitions.dimension_definitions import angle

# Description
## Most of cases might be represented in 2-dimensional space with two orthogonal axis - vertical Y and horizontal X. Any vector in this space (velocity, force etc) can be easily
## transformed to sum of two orthogonal vectors - projections of this vector to axis with help of angle between the vector and the axis.
## The best way to choose axis is to make gravitation vector colinear with Y axis, this zeroes X-projection of such vector.

vector_length  = symbols('vector_length')
vector_angle = symbols('vector_angle')
projection = symbols('projection')

law = Eq(projection, vector_length * cos(vector_angle))

def print():
    return pretty(projection, use_unicode=False)

#@validate_input(angle_ = angle)
#@validate_output()
def get_projection(vector_: Quantity, angle_: Quantity) -> Quantity:
    result_projection_expr = solve(law, projection, dict=True)[0][projection]
    result_expr = result_projection_expr.subs({vector_length: vector_, vector_angle: angle_})
    return expr_to_quantity(result_expr, 'projection_to_axis')

