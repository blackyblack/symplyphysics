from sympy import (Eq, Derivative, dsolve)
from symplyphysics import (units, Quantity, Function, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import mass_flow_rate as mass_flow_rate_law
from symplyphysics.definitions import density_from_mass_volume as density_law

# Description
## The rate of change in the volume of matter. For example, the outflow of a substance
## from a certain volume, the flow in a pipe section, the combustion of fuel.

# Definition: v = dV / dt
# Where:
## v is volume flow rate
## V is volume of matter (function V(t))
## t is time

time = Symbol("time", units.time)
volume = Function("volume", units.volume)
volume_flow_rate = Function("volume_flow_rate", units.volume / units.time)

definition = Eq(volume_flow_rate(time), Derivative(volume(time), time))

definition_units_SI = units.meters ** 3 / units.second


def print_law() -> str:
    return print_expression(definition)


@validate_input(volume_start_=volume, volume_end_=volume, time_=time)
@validate_output(volume_flow_rate)
def calculate_volume_flow_rate(volume_start_: Quantity, volume_end_: Quantity,
    time_: Quantity) -> Quantity:
    volume_function_ = time * (volume_end_ - volume_start_) / time_
    applied_definition = definition.subs(volume(time), volume_function_)
    # calculate volume flow rate
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
