from symplyphysics import (
    symbols, Eq, Function, simplify, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.definitions import power_is_energy_derivative as power_derivative

# Description
# Power directly proportional to energy (work) and inversely proportional to time
# P = Q / t
# where:
# Q  - some energy (work)
# t -  energy action time
# TODO: this law should be derived from definition power_is_energy_derivative
power, energy, time = symbols('power energy time')
law = Eq(power, energy / time)
constant_1, constant_2 = symbols('constant_1 constant_2')
linear_energy = constant_1 * power_derivative.time +constant_2
definition_applied = power_derivative.definition.rhs.subs(power_derivative.energy_function(time), linear_energy)
linear_power = definition_applied.doit()
definition_energy = linear_energy.subs({constant_1: power, constant_2: 0, power_derivative.time: time})
dsolved = definition_energy.doit()
result = Eq(energy, dsolved)
applied_law = solve(result, power, dict=True)[0][power]
# Check : derived power is same as declared
difference = simplify(applied_law - law.rhs)
assert(difference == 0)

def print():
    return pretty(law, use_unicode=False)

@validate_input(energy_=units.energy, time_=units.time)
@validate_output(units.power)
def calculate_power(energy_: Quantity, time_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({energy: energy_, time: time_})
    return expr_to_quantity(result_expr, 'power')