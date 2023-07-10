from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The instantaneous electrical current, or simply the electrical current, is the time derivative of the charge that flows.

# Definition: I = dQ/dt
# Where:
## I is current through the conductor
## Q is the electrical charge transferred

time = Symbol("time", units.time)
current = Function("current", units.current)
charge = Function("charge", units.charge)

definition = Eq(current(time), Derivative(charge(time), time))

definition_units_SI = units.ampere


def print_law() -> str:
    return print_expression(definition)


@validate_input(charge_start_=charge, charge_end_=charge, time_=time)
@validate_output(current)
def calculate_current(charge_start_: Quantity, charge_end_: Quantity, time_: Quantity) -> Quantity:
    charge_function_ = time * (charge_end_ - charge_start_) / time_
    applied_definition = definition.subs(charge(time), charge_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
