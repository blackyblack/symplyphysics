from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Current flowing through the resistor is proportional to applied voltage and reversly proportional to impedance of that resistor
## Ohm's law: I = U / R
## Where I is the current flowing through element
## U - is voltage drop on this element
## R is impedance of this element

current = Symbol("current", units.current)
voltage = Symbol("voltage", units.voltage)
resistance = Symbol("resistance", units.impedance)

law = Eq(current, voltage / resistance)

def print(expr: Expr) -> str:
    symbols = [current, voltage, resistance]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(voltage_=voltage, resistance_=resistance)
@validate_output_symbol(current)
def calculate_current(voltage_: Quantity, resistance_: Quantity) -> Quantity:
    result_current_expr = solve(law, current, dict=True)[0][current]
    result_expr = result_current_expr.subs({voltage: voltage_, resistance: resistance_})
    return expr_to_quantity(result_expr)
