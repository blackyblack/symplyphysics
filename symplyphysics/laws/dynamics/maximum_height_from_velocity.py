from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
# Description
## The maximum height to which a body thrown vertically upwards will rise depends on the initial velocity
## Law: h = (v**2)/(2*g)
## v is initial velocity,
## g is acceleration of free fall.

maximum_height = Symbol("maximum_height", units.length)
initial_velocity = Symbol("initial_velocity", units.velocity)

law = Eq(maximum_height, initial_velocity**2 / (2 * units.acceleration_due_to_gravity))


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_velocity_=initial_velocity)
@validate_output(maximum_height)
def calculate_maximum_height(initial_velocity_: Quantity) -> Quantity:
    result_maximum_height = solve(law, maximum_height, dict=True)[0][maximum_height]
    result_expr = result_maximum_height.subs({initial_velocity: initial_velocity_})
    return Quantity(result_expr)
