from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

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

def print(expr: Expr) -> str:
    symbols = [power, current, voltage]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(current_=current, voltage_=voltage)
@validate_output_symbol(power)
def calculate_power(current_: Quantity, voltage_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({current: current_, voltage: voltage_})
    return expr_to_quantity(result_expr)
