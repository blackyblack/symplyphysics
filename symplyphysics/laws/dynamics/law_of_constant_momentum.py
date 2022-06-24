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

