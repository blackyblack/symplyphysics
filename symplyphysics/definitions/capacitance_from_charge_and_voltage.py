from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## The electrical capacitance of a capacitor is
## charge accumulated in this capacitor divided to voltage on it

## Definition: C = Q / U
## Where:
## C is capacitance
## Q is charge
## U is voltage

capacitance, charge, voltage = symbols('capacitance charge voltage')
definition = Eq(capacitance, charge / voltage)

definition_dimension_SI = units.farad

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(charge_=units.charge, voltage_=units.voltage)
@validate_output(units.capacitance)
def calculate_capacitance(charge_: Quantity, voltage_: Quantity) -> Quantity:
    solved = solve(definition, capacitance, dict=True)[0][capacitance]
    result_expr = solved.subs({
        charge: charge_,
        voltage: voltage_})
    return expr_to_quantity(result_expr, 'capacitance')

@validate_input(capacitance_=units.capacitance, charge_=units.charge)
@validate_output(units.voltage)
def calculate_voltage(capacitance_: Quantity, charge_: Quantity) -> Quantity:
    solved = solve(definition, voltage, dict=True)[0][voltage]
    result_expr = solved.subs({
        capacitance: capacitance_,
        charge: charge_})
    return expr_to_quantity(result_expr, 'voltage')

@validate_input(capacitance_=units.capacitance, voltage_=units.voltage)
@validate_output(units.capacitance)
def calculate_charge(capacitance_: Quantity, voltage_: Quantity) -> Quantity:
    solved = solve(definition, charge, dict=True)[0][charge]
    result_expr = solved.subs({
        capacitance: capacitance_,
        voltage: voltage_})
    return expr_to_quantity(result_expr, 'accumulated_charge')