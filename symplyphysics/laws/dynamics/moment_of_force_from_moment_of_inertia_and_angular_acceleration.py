from sympy import Eq, solve
from symplyphysics import (angle_type, units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The moment of force (the moment of force relative to a point)
## is a vector physical quantity that characterizes the effect of force on a mechanical object,
## which can cause its rotational movement.

# Law: M = I * epsilon
# Where:
## epsilon is angular acceleration
## I - body moment of inertia
## M - moment of force

angular_acceleration = Symbol("angular_acceleration", angle_type / (units.time**2))
moment_of_inertia = Symbol("moment_of_inertia", units.mass * units.area)
moment_of_force = Symbol("moment_of_force", units.force * units.length * angle_type)

law = Eq(moment_of_force, moment_of_inertia * angular_acceleration)


def print_law() -> str:
    return print_expression(law)


@validate_input(moment_of_inertia_=moment_of_inertia, angular_acceleration_=angular_acceleration)
@validate_output(moment_of_force)
def calculate_moment_of_force(moment_of_inertia_: Quantity,
    angular_acceleration_: Quantity) -> Quantity:
    solved = solve(law, moment_of_force, dict=True)[0][moment_of_force]
    result_expr = solved.subs({
        moment_of_inertia: moment_of_inertia_,
        angular_acceleration: angular_acceleration_
    })
    return Quantity(result_expr)
