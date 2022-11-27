from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Dissipated heat power is proportional to current square and resistance
## Joule-Lenz law: P = I**2 * R
## P - is dissipated heat power of the element
## I is current flowing through this element
## R is impedance of this element

heat_power, current, resistance = symbols('heat_power current resistance')
law = Eq(heat_power, current**2 * resistance)

def print():
    return pretty(law, use_unicode=False)

@validate_input(current_=units.current, resistance_=units.impedance)
@validate_output(units.power)
def calculate_heat_power(current_: Quantity, resistance_: Quantity) -> Quantity:
    result_power_expr = solve(law, heat_power, dict=True)[0][heat_power]
    result_expr = result_power_expr.subs({current: current_, resistance: resistance_})
    return expr_to_quantity(result_expr, 'heat_power')
