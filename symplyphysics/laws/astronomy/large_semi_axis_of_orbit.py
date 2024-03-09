from sympy import (Eq, solve,)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)
from sympy.physics.units import gravitational_constant

# Description
## Let the body move in an elliptical orbit. Then its large semi-axis depends on the mass of the body around which it
## rotates and the orbital velocity.

## Law is: a = G * M / v^2, where
## a - large semi-axis of the orbit,
## G - gravitational constant,
## M - mass of the body around which the movement takes place,
## v - orbital velocity.

large_semi_axis = Symbol("large_semi_axis", units.length)

orbital_velocity = Symbol("orbital_velocity", units.velocity)
mass = Symbol("mass", units.mass)

law = Eq(large_semi_axis, gravitational_constant * mass / orbital_velocity**2)


def print_law() -> str:
    return print_expression(law)


@validate_input(orbital_velocity_=orbital_velocity, mass_=mass)
@validate_output(large_semi_axis)
def calculate_large_semi_axis(orbital_velocity_: Quantity, mass_: Quantity) -> Quantity:
    result_expr = solve(law, large_semi_axis, dict=True)[0][large_semi_axis]
    result_expr = result_expr.subs({
        orbital_velocity: orbital_velocity_,
        mass: mass_,
    })
    return Quantity(result_expr)
