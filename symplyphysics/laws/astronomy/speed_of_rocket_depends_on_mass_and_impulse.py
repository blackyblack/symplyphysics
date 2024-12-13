from sympy import (Eq, solve, log)
from symplyphysics import (
    clone_as_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
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

# Links: Wikipedia <https://en.wikipedia.org/wiki/Tsiolkovsky_rocket_equation>

speed = Symbol("speed", units.velocity)

specific_impulse = Symbol("specific_impulse", units.velocity)
initial_mass = clone_as_symbol(symbols.mass, display_symbol="m_0", display_latex="m_0")
final_mass = clone_as_symbol(symbols.mass, display_symbol="m_1", display_latex="m_1")

law = Eq(speed, specific_impulse * log(initial_mass / final_mass))


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
