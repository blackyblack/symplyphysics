from sympy import Expr
from symplyphysics import (
     Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Newton's second law: a = F / m

force = Symbol("force", units.force)
mass = Symbol("mass", units.mass)
acceleration = Symbol("acceleration", units.acceleration)

law = Eq(acceleration, force / mass)

def print(expr: Expr) -> str:
    symbols = [force, mass, acceleration]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(mass_=mass, acceleration_=acceleration)
@validate_output_symbol(force)
def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    result_force_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_force_expr.subs({mass: mass_, acceleration: acceleration_})
    return expr_to_quantity(result_expr)
