from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

# Description
# Amount of energy for body heat Q = C * m * (t2 - t1)
# where:
# C - specific heat capacity
# m - body mass
# t1 - initial temperature
# t2 - final temperature

# Note: specific heat capacity C is parameter indicates how much energy in joules
# must be expended to heat a unit mass of a substance by one unit temperature.

# Note: the resultant energy is the energy consumed for changing body temperature from value t1 to value t2
# Positive value Q is for absorbed energy (heating), negative for released energy (cooling)

amount_energy = Symbol("amount_energy", units.energy)
specific_heat_capacity = Symbol("specific_heat_capacity",
    units.energy / (units.mass * units.temperature))
temperature_origin = Symbol("temperature_origin", units.temperature)
temperature_end = Symbol("temperature_end", units.temperature)

law = Eq(amount_energy,
    specific_heat_capacity * symbols.basic.mass * (temperature_end - temperature_origin))


def print_law() -> str:
    return print_expression(law)


@validate_input(specific_heat_capacity_=specific_heat_capacity,
    body_mass_=symbols.basic.mass,
    temperature_end_=temperature_end,
    temperature_origin_=temperature_origin)
@validate_output(amount_energy)
def calculate_amount_energy(specific_heat_capacity_: Quantity, body_mass_: Quantity,
    temperature_end_: Quantity, temperature_origin_: Quantity) -> Quantity:

    result_amount_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_amount_energy_expr.subs({
        specific_heat_capacity: specific_heat_capacity_,
        symbols.basic.mass: body_mass_,
        temperature_end: temperature_end_,
        temperature_origin: temperature_origin_
    })
    return Quantity(result_expr)
