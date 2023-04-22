from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## The electrical capacitance of a capacitor is
## charge accumulated in this capacitor divided to voltage on it

# Definition: C = Q / U
# Where:
## C is capacitance
## Q is charge
## U is voltage

capacitance = Symbol("capacitance", units.capacitance)
charge = Symbol("charge", units.charge)
voltage = Symbol("voltage", units.voltage)

definition = Eq(capacitance, charge / voltage)

definition_units_SI = units.farad


def print() -> str:
    return print_expression(definition)


@validate_input_symbols(charge_=charge, voltage_=voltage)
@validate_output_symbol(capacitance)
def calculate_capacitance(charge_: Quantity, voltage_: Quantity) -> Quantity:
    solved = solve(definition, capacitance, dict=True)[0][capacitance]
    result_expr = solved.subs({charge: charge_, voltage: voltage_})
    return expr_to_quantity(result_expr)
