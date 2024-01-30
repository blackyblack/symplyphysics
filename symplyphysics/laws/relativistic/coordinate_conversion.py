from sympy import (Eq, solve, sqrt)
from sympy.physics.units import speed_of_light
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)

# Description
## Let the frame of reference move relative to the other frame of reference at a constant speed along the X
## axis, and the origin of the spatial coordinates coincide at the initial moment of time in both systems.
## Then there are simple transformations that can be used to get the x coordinate in one frame of reference,
## knowing the x coordinate in another frame of reference.

## Law is: x_2 = (x_1 - v * t) / (1 - (v / c)^2), where
## x_2 - another coordinate (in the second frame of reference),
## x_1 - coordinate (in the first frame of reference),
## v - velocity at which the second frame of reference moves relative to the first,
## a - time (in the first frame of reference).

another_coordinate = Symbol("another_coordinate", units.length)

coordinate = Symbol("coordinate", units.length)
velocity = Symbol("velocity", units.velocity)
time = Symbol("time", units.time)

law = Eq(another_coordinate, (coordinate - velocity * time) / (1 - (velocity / speed_of_light)**2))

def print_law() -> str:
    return print_expression(law)


@validate_input(coordinate_=coordinate, velocity_=velocity, time_=time)
@validate_output(another_coordinate)
def calculate_another_coordinate(coordinate_: Quantity, velocity_: Quantity, time_: Quantity) -> Quantity:
    result_another_coordinate_expr = solve(law, another_coordinate, dict=True)[0][another_coordinate]
    result_expr = result_another_coordinate_expr.subs({coordinate: coordinate_, velocity: velocity_, time: time_})
    return Quantity(result_expr)
