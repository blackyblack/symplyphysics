from sympy import (Eq, solve, log)
from symplyphysics import (
    clone_symbol,
    symbols,
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

## Law is: V = I * ln(M1 / M2), where
## V - final speed of the rocket,
## I- specific impulse of a rocket engine,
## M1 - initial mass of the rocket,
## M2 - final mass of the rocket.

speed = Symbol("speed", units.velocity)

specific_impulse = Symbol("specific_impulse", units.velocity)
initial_mass = symbols.mass
final_mass = symbols.mass

law = Eq(speed, specific_impulse * log(initial_mass / final_mass))


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
