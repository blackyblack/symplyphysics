from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Self-induction voltage (or EMF - electromotive force) can be expressed in form of modified Faraday's law.
## In a circuit, it is much easier to measure current than it is to measure magnetic flux, so the induced voltage
## can be found if the inductance and frequency of the current are known.

## Self-induction voltage law: E = -L * dI/dt, where
## E is the self-induction voltage
## L is the inductance of the inductor
## I(t) is the current through the inductor

time = Symbol("time", units.time)
self_induction_voltage = Function("self_induction_voltage", units.voltage)
current = Function("current", units.current)
inductance = Symbol("inductance", units.inductance)

definition = Eq(self_induction_voltage(time), -1 * inductance * Derivative(current(time), time))

definition_units_SI = units.volt


def print_law() -> str:
    return print_expression(definition)


@validate_input(inductance_=inductance, current_start_=current, current_end_=current, time_=time)
@validate_output(self_induction_voltage)
def calculate_voltage(inductance_: Quantity, current_start_: Quantity, current_end_: Quantity,
    time_: Quantity) -> Quantity:
    current_function_ = time * (current_end_ - current_start_) / time_
    applied_definition = definition.subs({
        current(time): current_function_,
        inductance: inductance_
    })
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
