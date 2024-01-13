from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The rate of change in the mass of matter. For example, the outflow of a substance
## from a certain volume, the flow in a pipe section, the combustion of fuel.

# Definition: mu = dm / dt
# Where:
## mu is mass flow rate
## t is time

time = Symbol("time", units.time)
mass = Function("mass", units.mass)
mass_flow_rate = Function("mass_flow_rate", units.mass / units.time)

definition = Eq(mass_flow_rate(time), Derivative(mass(time), time))

definition_units_SI = units.kilogram / units.second


def print_law() -> str:
    return print_expression(definition)


@validate_input(mass_start_=mass, mass_end_=mass, time_=time)
@validate_output(mass_flow_rate)
def calculate_mass_flow_rate(mass_start_: Quantity, mass_end_: Quantity,
    time_: Quantity) -> Quantity:
    mass_function_ = time * (mass_end_ - mass_start_) / time_
    applied_definition = definition.subs(mass(time), mass_function_)
    # calculate mass flow rate
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
