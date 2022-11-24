from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
# Power of current is proportional current and voltage
# P = I * U
# where:
# I -the current flowing through element
# U - is voltage  on this element
operate_power, operate_current, operate_voltage = symbols('operate_power operate_current operate_voltage')
law = Eq(operate_power, operate_current * operate_voltage)

def print():
    return pretty(law, use_unicode=False)

@validate_input(operate_current_=units.current, operate_voltage_=units.voltage)
@validate_output(units.power)
def calculate_power(operate_current_: Quantity, operate_voltage_: Quantity) -> Quantity:
    result_power_expr = solve(law, operate_power, dict=True)[0][operate_power]
    result_expr = result_power_expr.subs({operate_current: operate_current_, operate_voltage: operate_voltage_})
    return expr_to_quantity(result_expr, 'operate_power')