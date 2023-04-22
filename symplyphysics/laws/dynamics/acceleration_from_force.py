from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol,
                           print_expression, validate_input_symbols,
                           validate_output_symbol)

# Description
## Newton's second law: a = F / m

force = Symbol("force", units.force)
mass = Symbol("mass", units.mass)
acceleration = Symbol("acceleration", units.acceleration)

law = Eq(acceleration, force / mass)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(mass_=mass, acceleration_=acceleration)
@validate_output_symbol(force)
def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    result_force_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_force_expr.subs({
        mass: mass_,
        acceleration: acceleration_
    })
    return expr_to_quantity(result_expr)
