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
volume = Function("volume", units.mass)
volume_flow_rate = Function("volume_flow_rate", units.mass / units.time)

definition = Eq(volume_flow_rate(time), Derivative(volume(time), time))

definition_units_SI = units.meters ** 3 / units.second

# If the mass and volume are constant over time,
# then a standard formula can be obtained: pho = m / V
mass_rate_equation = mass_flow_rate_law.definition.subs({
    mass_flow_rate_law.time: time,
    mass_flow_rate_law.mass_flow_rate(mass_flow_rate_law.time): 0
})
mass_equation = dsolve(mass_rate_equation, mass_flow_rate_law.mass(time)).subs({
    "C1": density_law.mass
})

volume_rate_equation = definition.subs({
    volume_flow_rate(time): 0
})
volume_equation = dsolve(volume_rate_equation, volume(time)).subs({
    "C1": density_law.volume
})

density_equation = density_law.definition.subs({
    density_law.mass: mass_equation.rhs,
    density_law.volume: volume_equation.rhs
})
assert expr_equals(density_law.definition.rhs, density_equation.rhs)


def print_law() -> str:
    return print_expression(definition)


@validate_input(volume_start_=volume, volume_end_=volume, time_=time)
@validate_output(volume_flow_rate)
def calculate_mass_flow_rate(volume_start_: Quantity, volume_end_: Quantity,
    time_: Quantity) -> Quantity:
    mass_function_ = time * (volume_end_ - volume_start_) / time_
    applied_definition = definition.subs(volume(time), mass_function_)
    # calculate volume flow rate
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
