from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Ohm's law - I = U / R
## Current flowing through the resistor is proportional to applied voltage and reversly proportional to impedance of that resistor

current, voltage, resistance = symbols('current voltage resistance')
law = Eq(current, voltage / resistance)

def print():
    return pretty(law, use_unicode=False)

@validate_input(voltage_=units.voltage, resistance_=units.impedance)
@validate_output(units.current)
def calculate_current(voltage_: Quantity, resistance_: Quantity) -> Quantity:
    result_current_expr = solve(law, current, dict=True)[0][current]
    result_expr = result_current_expr.subs({voltage: voltage_, resistance: resistance_})
    return expr_to_quantity(result_expr, 'current')
