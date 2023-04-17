from sympy import Expr
from symplyphysics import (
    Eq, pretty, units, solve, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Wavelength is the spatial period of a periodic wave — the distance over which the wave's shape repeats.
## It is the distance between consecutive corresponding points of the same phase on the wave.
## Law: λ = v * T, where
## λ is wavelength,
## v is wave propagation speed (phase speed),
## T is oscillation period.

#TODO derive this from velocity and period definitions

wavelength = Symbol("wavelength", units.length)
propagation_speed = Symbol("propagation_speed", units.velocity)
oscillation_period = Symbol("oscillation_period", units.time)

law = Eq(wavelength, propagation_speed * oscillation_period)

def print(expr: Expr) -> str:
    symbols = [wavelength, propagation_speed, oscillation_period]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(velocity_=propagation_speed, period_=oscillation_period)
@validate_output_symbol(wavelength)
def calculate_wavelength(velocity_: Quantity, period_: Quantity) -> Quantity:
    applied_definition = solve(law, wavelength, dict=True)[0][wavelength]
    result_expr = applied_definition.subs({propagation_speed: velocity_, oscillation_period: period_})
    return expr_to_quantity(result_expr)
