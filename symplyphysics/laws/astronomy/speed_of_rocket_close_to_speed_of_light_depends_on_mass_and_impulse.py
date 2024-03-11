from sympy import (Eq, solve)
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
## The Tsiolkovsky formula determines the speed that an aircraft develops under the influence of the thrust
## of a rocket engine, unchanged in direction, in the absence of all other forces.
## For a rocket flying at a speed close to the speed of light, the generalized Tsiolkovsky formula is valid,
## in which the speed of light is present.

## Law is: M2 / M1 = ((1 - (V / c)) / (1 + (V / c))) ^(c / (2 * I)), where
## V - final speed of the rocket,
## I- specific impulse of a rocket engine,
## M1 - initial mass of the rocket,
## M2 - final mass of the rocket,
## c - speed of light.

speed = Symbol("speed", units.velocity)

specific_impulse = Symbol("specific_impulse", units.velocity)
initial_mass = Symbol("initial_mass", units.mass)
final_mass = Symbol("final_mass", units.mass)

law = Eq(final_mass / initial_mass, ((1 - (speed / speed_of_light)) / (1 + (speed / speed_of_light)))**(speed_of_light / (2 * specific_impulse)))


def print_law() -> str:
    return print_expression(law)


@validate_input(specific_impulse_=specific_impulse,
    initial_mass_=initial_mass,
    final_mass_=final_mass)
@validate_output(speed)
def calculate_speed(specific_impulse_: Quantity, initial_mass_: Quantity,
    final_mass_: Quantity) -> Quantity:
    result_expr = solve(law, speed, dict=True)[0][speed]
    result_expr = result_expr.subs({
        specific_impulse: specific_impulse_,
        initial_mass: initial_mass_,
        final_mass: final_mass_
    })
    return Quantity(result_expr)
