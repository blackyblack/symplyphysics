from symplyphysics import (
    symbols, Eq, simplify, pretty, solve, Quantity, units,
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
# Where :
# constant_1, constant_2 are formal parameters for the linear energy function condition
# constant_1 - constant_power
# constant_2 - initial energy.It is equal 0.
energy_applied = constant_1 * power_derivative.time +constant_2
energy_applied = energy_applied.subs(constant_2, 0)
definition_applied = power_derivative.definition.rhs.subs(power_derivative.energy_function(power_derivative.time), energy_applied)
power_applied = definition_applied.doit()
power_applied_law = Eq(power, power_applied)
energy_applied_law = Eq(energy, energy_applied)
applied_law = solve([power_applied_law,energy_applied_law], (power,constant_1), dict=True)[0][constant_1]
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