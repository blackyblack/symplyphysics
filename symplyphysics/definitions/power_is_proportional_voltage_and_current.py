from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
# Power of current is proportional to current and voltage
# P = I * U
# where:
# I -the current flowing through element
# U - is voltage  on this element
power, current, voltage = symbols('power current voltage')
law = Eq(power, current * voltage)

def print():
    return pretty(law, use_unicode=False)

@validate_input(current_=units.current, voltage_=units.voltage)
@validate_output(units.power)
def calculate_power(current_: Quantity, voltage_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({current: current_, voltage: voltage_})
    return expr_to_quantity(result_expr, 'power')