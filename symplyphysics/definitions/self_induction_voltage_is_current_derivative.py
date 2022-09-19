from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Self-induction voltage definition: E = -L * dI/dt, where
## E is the self-induction voltage
## L is the inductance of the inductor
## I(t) is the current through the inductor

time = symbols('time')
voltage_function, current_function = symbols('voltage current', cls = Function)
inductance = symbols('inductance')
definition = Eq(voltage_function(time), -1 * inductance * Derivative(current_function(time), time))

definition_dimension_SI = units.volt

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(inductance_=units.inductance, current_start_=units.current, current_end_=units.current, time_=units.time)
@validate_output(units.voltage)
def calculate_voltage(inductance_: Quantity, current_start_: Quantity, current_end_: Quantity, time_: Quantity) -> Quantity:
    current_function_ = (current_end_ - current_start_) / time_
    applied_definition = definition.subs(current_function(time), current_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, 'voltage')
