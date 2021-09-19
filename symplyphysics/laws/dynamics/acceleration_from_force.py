from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Newton's second law: a = F / m

force, mass, acceleration = symbols('force mass acceleration')
law = Eq(acceleration, force / mass)

def print():
    return pretty(law, use_unicode=False)

@validate_input(mass_=units.mass, acceleration_=units.acceleration)
@validate_output(units.force)
def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    result_expr = solve(law.subs({mass: mass_, acceleration: acceleration_}))[0]
    return expr_to_quantity(result_expr, 'force')
