from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
#If there is no external force is applied to system of objects, the summary momentum of this system remains constant during and after any interactions between objects
#Scalar version for two colliding objects in 1-dimentional env. m1, V1, m2, V2 are masses and Velocities befor collision, m3, V3, m4 and V4 - after.
#Also applicable for reactive engine simulation

momentum_before, momentum_after = symbols('momentum_before momentum_after')
law = Eq(momentum_after, momentum_before)

def print():
    return pretty(law, use_unicode=False)

@validate_input(momentum_before = units.kilogram * units.meter / units.second, momentum_after = units.kilogram * units.meter / units.second)
@validate_output(units.momentum)

def calculate_momentum_after(momentum_before_: Quantity) -> Quantity:
    solved = solve(law, momentum_before, dict=True)[0][momentum_after]
    result_expr = solved.subs()
    return expr_to_quantity(result_expr, 'momentum after')
