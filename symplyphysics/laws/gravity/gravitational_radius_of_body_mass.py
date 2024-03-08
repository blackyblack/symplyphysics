from sympy import (Eq, solve,)
from sympy.physics.units import gravitational_constant, speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The gravitational radius is a characteristic radius defined for any physical body with mass. This is the radius of the sphere
## on which the event horizon created by this mass would be located (from the point of view of general theory of relativity) if it were distributed spherically
## symmetrically, would be stationary (in particular, it would not rotate, but radial movements are permissible) and would lie entirely
## inside this sphere.

## Law is: R = 2 * G * M / c^2, where
## R - gravitational radius of the body,
## G - gravitational constant,
## M - mass of the body,
## c - speed of light.

gravitational_radius = Symbol("gravitational_radius", units.length)

body_mass = Symbol("body_mass", units.mass)

law = Eq(gravitational_radius, 2 * gravitational_constant * body_mass / speed_of_light**2)


def print_law() -> str:
    return print_expression(law)


@validate_input(body_mass_=body_mass)
@validate_output(gravitational_radius)
def calculate_radius(body_mass_: Quantity) -> Quantity:
    result_expr = solve(law, gravitational_radius, dict=True)[0][gravitational_radius]
    result_expr = result_expr.subs({
        body_mass: body_mass_,
    })
    return Quantity(result_expr)
