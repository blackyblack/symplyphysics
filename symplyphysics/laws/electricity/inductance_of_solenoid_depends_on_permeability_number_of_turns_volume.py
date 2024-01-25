from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)
from sympy.physics.units import magnetic_constant

# Description
## A solenoid is a cylindrical coil consisting of a large number of turns of wire forming a helical line.
## Inductance of solenoid depends on core material, number of turns per unit length, and volume of core.

## Law is: L = mu * mu0 * n^2 * V, where
## L - inductance,
## mu - relative permeability of the core inside of a solenoid,
## mu0 - magnetic constant,
## n - number of turns per unit length,
## V - volume of solenoid.

inductance = Symbol("inductance", units.inductance)

relative_permeability = Symbol("relative_permeability", dimensionless)
number_of_turns = Symbol("number_of_turns", 1 / units.length)
volume = Symbol("volume", units.volume)

law = Eq(inductance, relative_permeability * magnetic_constant * number_of_turns**2 * volume)


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permeability_=relative_permeability, number_of_turns_=number_of_turns, volume_=volume)
@validate_output(inductance)
def calculate_inductance(relative_permeability_: float, number_of_turns_: Quantity, volume_: Quantity) -> Quantity:
    result_inductance_expr = solve(law, inductance, dict=True)[0][inductance]
    result_expr = result_inductance_expr.subs({relative_permeability: relative_permeability_, number_of_turns: number_of_turns_, volume: volume_})
    return Quantity(result_expr)
