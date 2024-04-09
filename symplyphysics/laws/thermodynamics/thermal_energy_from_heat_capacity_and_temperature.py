from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
# Amount of energy for body heat Q = C * (t2 - t1)
# where:
# C - heat capacity
# t1 - initial temperature
# t2 - final temperature

# Note: heat capacity C is parameter indicates how much energy in joules
# must be expended to heat the substance by one unit temperature.

# Note: the resultant energy is the energy consumed for changing body temperature from value t1 to value t2
# Positive value Q is for absorbed energy (heating), negative for released energy (cooling)

# Note: most of the time, intensive heat capacity is known. Use [molar](./heat_capacity_via_molar_heat_capacity.py)
# and [specific](heat_capacity_via_specific_heat_capacity.py) laws for this.

amount_energy = Symbol("amount_energy", units.energy)
heat_capacity = Symbol("heat_capacity", units.energy / units.temperature)
temperature_origin = Symbol("temperature_origin", units.temperature)
temperature_end = Symbol("temperature_end", units.temperature)

law = Eq(amount_energy, heat_capacity * (temperature_end - temperature_origin))


def print_law() -> str:
    return print_expression(law)


@validate_input(heat_capacity_=heat_capacity,
    temperature_end_=temperature_end,
    temperature_origin_=temperature_origin)
@validate_output(amount_energy)
def calculate_amount_energy(heat_capacity_: Quantity,
    temperature_end_: Quantity, temperature_origin_: Quantity) -> Quantity:

    result_amount_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_amount_energy_expr.subs({
        heat_capacity: heat_capacity_,
        temperature_end: temperature_end_,
        temperature_origin: temperature_origin_
    })
    return Quantity(result_expr)
