from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity, sqrt
)
from sympy.physics.units import speed_of_light, magnetic_constant, electric_constant
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol

from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Speed of light in vacuum is fundamental but still might be calculated from other fundamentals.

# Law: c = 1/sqrt(e0 * u0), where
## C is speed of light in vacuum,
## e0 is electric constant or vacuum permittivity,
## u0 is magnetic constant or vacuum permeability.


law = Eq(speed_of_light, 1 / sqrt(magnetic_constant * electric_constant))

def print(expr: Expr) -> str:
    symbols = [speed_of_light, magnetic_constant, electric_constant]
    return pretty(to_printable(expr, symbols), use_unicode=False)

assert speed_of_light == 1 / sqrt(magnetic_constant * electric_constant)

