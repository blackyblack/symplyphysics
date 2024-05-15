from sympy import (Eq, solve, sqrt)
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Let the frame of reference move relative to the other frame of reference at a constant speed along the X
## axis, and the origin of the spatial coordinates coincide at the initial moment of time in both systems.
## Then there are simple transformations that can be used to get the x coordinate in one frame of reference,
## knowing the x coordinate in another frame of reference.

## Law is: x_2 = (x_1 - v * t) / sqrt(1 - (v / c)^2), where
## x_2 - another coordinate (in the second frame of reference),
## x_1 - coordinate (in the first frame of reference),
## v - velocity at which the second frame of reference moves relative to the first,
## t - time in the first frame of reference.

coordinate_second_frame = Symbol("coordinate_second_frame", units.length)

coordinate_first_frame = Symbol("coordinate_first_frame", units.length)
velocity = Symbol("velocity", units.velocity)
time_first_frame = Symbol("time_first_frame", units.time)

law = Eq(coordinate_second_frame, (coordinate_first_frame - velocity * time_first_frame) / sqrt(
    (1 - (velocity / speed_of_light)**2)))


def print_law() -> str:
    return print_expression(law)


@validate_input(coordinate_first_frame_=coordinate_first_frame,
    velocity_=velocity,
    time_first_frame_=time_first_frame)
@validate_output(coordinate_second_frame)
def calculate_coordinate_second_frame(coordinate_first_frame_: Quantity, velocity_: Quantity,
    time_first_frame_: Quantity) -> Quantity:
    result_coordinate_first_frame_second_frame_expr = solve(law, coordinate_second_frame,
        dict=True)[0][coordinate_second_frame]
    result_expr = result_coordinate_first_frame_second_frame_expr.subs({
        coordinate_first_frame: coordinate_first_frame_,
        velocity: velocity_,
        time_first_frame: time_first_frame_
    })
    return Quantity(result_expr)
