from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Inductor can accumlate energy in the magnetic field inside it.
## Law: Q = L * I^2 / 2
## Q is accumulated energy
## L is inductance of inductor
## I is current flowing through the inductor

accumulated_energy = Symbol("accumulated_energy", units.energy)
inductance = Symbol("inductance", units.inductance)
current = Symbol("current", units.current)

law = Eq(accumulated_energy, inductance * current**2 / 2)

def print(expr: Expr) -> str:
    symbols = [accumulated_energy, inductance, current]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(inductance_=inductance, current_=current)
@validate_output_symbol(accumulated_energy)
def calculate_accumulated_energy(inductance_: Quantity, current_: Quantity) -> Quantity:
    result_energy_expr = solve(law, accumulated_energy, dict=True)[0][accumulated_energy]
    result_expr = result_energy_expr.subs({inductance: inductance_, current: current_})
    return expr_to_quantity(result_expr)
