from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Dimensionless, Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Volume number density is the number of specified objects per unit volume.

# Definition: N = n / V
# Where:
## n is the total number of objects
## V is volume

number_density = Symbol("number_density", 1 / units.volume)
objects = Symbol("objects", Dimensionless)
volume = Symbol("volume", units.volume)

definition = Eq(number_density, objects / volume)

definition_units_SI = 1 / units.meter**3

def print(expr: Expr) -> str:
    symbols = [number_density, objects, volume]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(objects_=objects, volume_=volume)
@validate_output_symbol(number_density)
def calculate_number_density(objects_: int, volume_: Quantity) -> Quantity:
    solved = solve(definition, number_density, dict=True)[0][number_density]
    result_expr = solved.subs({
        objects: objects_,
        volume: volume_})
    return expr_to_quantity(result_expr)
