from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

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
specific_heat_capacity = Symbol("specific_heat_capacity", units.energy / (units.mass * units.temperature))
body_mass = Symbol("body_mass", units.mass)
temperature_origin = Symbol("temperature_origin", units.temperature)
temperature_end = Symbol("temperature_end", units.temperature)

law = Eq(amount_energy, specific_heat_capacity * body_mass * (temperature_end - temperature_origin))

def print(expr: Expr) -> str:
    symbols = [amount_energy, specific_heat_capacity, body_mass, temperature_origin, temperature_end]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(
    specific_heat_capacity_=specific_heat_capacity,
    body_mass_=body_mass,
    temperature_end_=temperature_end,
    temperature_origin_=temperature_origin)
@validate_output_symbol(amount_energy)
def calculate_amount_energy(
    specific_heat_capacity_: Quantity,
    body_mass_: Quantity,
    temperature_end_: Quantity,
    temperature_origin_: Quantity) -> Quantity:

    result_amount_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_amount_energy_expr.subs({
        specific_heat_capacity: specific_heat_capacity_,
        body_mass: body_mass_,
        temperature_end: temperature_end_,
        temperature_origin: temperature_origin_})
    return expr_to_quantity(result_expr)
