from sympy import (Eq, solve, atan)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## A gravitational maneuver is a purposeful change in the trajectory and flight speed of a spacecraft under the influence
## of the gravitational fields of celestial bodies.
## The maximum angle of rotation of the rocket during a gravitational maneuver depends on the first cosmic velocity of
## the planet around which the maneuver is performed, and the speed of the rocket relative to this planet.
## Law describes the closest possible maneuver, when aiming range is at planet radius.

## Law is: phi = arctg((v1 / v2)^2), where
## phi - maximum angle of rotation during a gravitational maneuver (angle at which the velocity vector of the rocket rotates),
## v1 - first cosmic velocity for a given planet,
## v2 - rocket's velocity relative to the planet.

# Conditions:
# - law describes the closest possible maneuver, when aiming range is at planet radius.

maximum_angle = Symbol("maximum_angle", angle_type)

first_cosmic_velocity_planet= Symbol("first_cosmic_velocity_planet", units.velocity)
rocket_speed = Symbol("rocket_speed", units.velocity)

law = Eq(maximum_angle, atan((first_cosmic_velocity_planet/ rocket_speed)**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(first_cosmic_velocity_planet_=first_cosmic_velocity_planet, rocket_speed_=rocket_speed)
@validate_output(maximum_angle)
def calculate_maximum_angle(first_cosmic_velocity_planet_: Quantity, rocket_speed_: Quantity) -> Quantity:
    result_expr = solve(law, maximum_angle, dict=True)[0][maximum_angle]
    result_expr = result_expr.subs({
        first_cosmic_velocity_planet: first_cosmic_velocity_planet_,
        rocket_speed: rocket_speed_,
    })
    return Quantity(result_expr)
