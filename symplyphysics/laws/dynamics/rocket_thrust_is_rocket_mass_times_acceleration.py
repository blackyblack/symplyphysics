from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## ...

# Law: R * v_rel = M * a
## R - rate of fuel consumption
## v_rel - velocity of racket relative to products
## M - racket mass
## a - racket acceleration

# Note: R*v_rel is called thrust of rocket engine

fuel_consumption_rate = Symbol("fuel_consumption_rate", units.mass / units.time)
relative_velocity = Symbol("relative_velocity", units.velocity)
racket_mass = Symbol("racket_mass", units.mass)
racket_acceleration = Symbol("racket_acceleration", units.acceleration)

law = Eq(fuel_consumption_rate * relative_velocity, racket_mass * racket_acceleration)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    fuel_consumption_rate_=fuel_consumption_rate,
    racket_mass_=racket_mass,
    racket_acceleration_=racket_acceleration,
)
@validate_output(relative_velocity)
def calculate_relative_velocity(
    fuel_consumption_rate_: Quantity,
    racket_mass_: Quantity,
    racket_acceleration_: Quantity,
) -> Quantity:
    result = solve(law, relative_velocity)[0].subs({
        fuel_consumption_rate: fuel_consumption_rate,
        racket_mass: racket_mass_,
        racket_acceleration: racket_acceleration_,
    })
    return Quantity(result)
