from sympy import Eq, solve
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## Let the body move in an elliptical orbit. Then its large semi-axis depends on the mass of the body around which it
## rotates and the orbital velocity.

## Law is: a = G * M / v^2, where
## a - large half-axis of the orbit,
## G - gravitational constant,
## M - mass of the body around which the movement takes place,
## v - orbital velocity.

# TODO: find link

large_half_axis_length = Symbol("large_half_axis_length", units.length)
orbital_velocity = Symbol("orbital_velocity", units.velocity)
planet_mass = symbols.mass

law = Eq(large_half_axis_length, gravitational_constant * planet_mass / orbital_velocity**2)


@validate_input(orbital_velocity_=orbital_velocity, planet_mass_=planet_mass)
@validate_output(large_half_axis_length)
def calculate_large_half_axis_length(orbital_velocity_: Quantity,
    planet_mass_: Quantity) -> Quantity:
    result_expr = solve(law, large_half_axis_length, dict=True)[0][large_half_axis_length]
    result_expr = result_expr.subs({
        orbital_velocity: orbital_velocity_,
        planet_mass: planet_mass_,
    })
    return Quantity(result_expr)
