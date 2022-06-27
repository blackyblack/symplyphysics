from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## If there is no external force is applied to system of objects, the summary momentum of this system remains constant during and after any interactions between objects

mass1, mass2, velocity1, velocity2 = symbols('mass1 mass2 velocity1 velocity2')
law = Eq(mass2 * velocity2, mass1 * velocity1)

def print():
    return pretty(law, use_unicode=False)

@validate_input(mass1_=units.mass, mass2_=units.mass, valocity1_=units.velocity, velocity2_=units.velocity)
@validate_output(units.mass)
@validate_output(units.velocity)

def calculate_velocity(mass1_: Quantity, velocity1_: Quantity, mass2_ : Quantity) -> Quantity:
    result_velocity_expr = solve(law, mass1_, velocity1_, mass2_, dict=True)[0]
    result_expr = result_velocity_expr.subs({mass1: mass1_, velocity1: velocity1_, mass2: mass2_})
    return expr_to_quantity(result_expr, 'Velocity2')

def calculate_mass(mass1_: Quantity, velocity1_: Quantity, velocity2_: Quantity) -> Quantity:
    result_mass_expr = solve(law, mass1, velocity1, velocity2, dict=True)[0]
    result_expr = result_mass_expr.subs({mass1: mass1_, velocity1: velocity1_, velocity2: velocity2_})
    return expr_to_quantity(result_expr, 'mass2')