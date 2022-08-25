from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Electrical current definition: I = dQ/dt
## Where I is current through the conductor, Q is the electrical charge transferred

time = symbols('time')
current_function, charge_function = symbols('current charge', cls = Function)
definition = Eq(current_function(time), Derivative(charge_function(time), time))

definition_dimension_SI = units.ampere

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(charge_start_=units.charge, charge_end_=units.charge, time_=units.time)
@validate_output(units.current)
def calculate_current(charge_start_: Quantity, charge_end_: Quantity, time_: Quantity) -> Quantity:
    charge_function_ = time * (charge_end_ - charge_start_) / time_
    applied_definition = definition.subs(charge_function(time), charge_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, 'current')

@validate_input(current_start_=units.current, current_end_=units.current, time_=units.time)
@validate_output(units.charge)
def calculate_charge(current_start_: Quantity, current_end_: Quantity, time_: Quantity) -> Quantity:
    current_function_ = 0.5 * (current_end_ + current_start_) * time_ * time
    applied_definition = definition.subs(current_function(time), current_function_)
    solved = applied_definition.doit()
    result_expr = solved.rhs
    return expr_to_quantity(result_expr, 'charge')
