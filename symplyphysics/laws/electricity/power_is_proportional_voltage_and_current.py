from sympy import (Eq, solve)
from symplyphysics import (
    units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol
)

# Description
# Power of current is proportional to current and voltage
# P = I * U
# where:
# I - the current flowing through element
# U - is voltage  on this element

power = Symbol("power", units.power)
current = Symbol("current", units.current)
voltage = Symbol("voltage", units.voltage)

law = Eq(power, current * voltage)

def print() -> str:
    return print_expression(law)

@validate_input_symbols(current_=current, voltage_=voltage)
@validate_output_symbol(power)
def calculate_power(current_: Quantity, voltage_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({current: current_, voltage: voltage_})
    return expr_to_quantity(result_expr)
