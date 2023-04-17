from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Ideal gas law: P * V = n * R * T
## Where:
## P is pressure,
## V is volume,
## n is number of moles,
## R is ideal gas constant,
## T is temperature

pressure = Symbol("pressure", units.pressure)
volume = Symbol("volume", units.volume)
mole_count = Symbol("mole_count", units.amount_of_substance)
temperature = Symbol("temperature", units.temperature)

law = Eq(pressure, mole_count * temperature * units.molar_gas_constant / volume)

def print(expr: Expr) -> str:
    symbols = [pressure, volume, mole_count, temperature]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(volume_=volume, temperature_=temperature, mole_count_=mole_count)
@validate_output_symbol(pressure)
def calculate_pressure(volume_: Quantity, temperature_: Quantity, mole_count_: Quantity) -> Quantity:
    solved = solve(law, pressure, dict=True)[0][pressure]
    result_expr = solved.subs({
        volume: volume_,
        temperature: temperature_,
        mole_count: mole_count_})
    return expr_to_quantity(result_expr)
