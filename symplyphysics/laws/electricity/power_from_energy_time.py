from symplyphysics import (
    symbols, Eq, simplify, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.definitions import power_is_energy_derivative as power_derivative

# Description
## Power directly proportional to energy (work) and inversely proportional to time
## P = Q / t
## Where:
## Q - some energy (work)
## t - energy action time

# Conditions:
## - Energy is linear function, or power is constant
## - Initial energy is 0

power, energy, time = symbols("power energy time")
law = Eq(power, energy / time)

# Derive the same law from definition of power as energy derivative

linear_function_coefficient, initial_energy_constant = symbols("linear_function_coefficient initial_energy_constant")

energy_linear_function = linear_function_coefficient * time + initial_energy_constant

power_definition_applied = power_derivative.definition.subs({
    power_derivative.energy_function(power_derivative.time): energy_linear_function,
    power_derivative.power_function(power_derivative.time): power})

power_applied_eq = power_definition_applied.doit()
energy_eq = Eq(energy, energy_linear_function)
derived_power = solve([power_applied_eq, energy_eq], (power, linear_function_coefficient), dict=True)[0][power]

# derived_power = (energy - initial_energy_constant)/time
# Assume initial_energy_constant = 0, according to law conditions
derived_power_without_initial_energy = derived_power.subs(initial_energy_constant, 0)

# Check that derived power is same as declared
assert(simplify(derived_power_without_initial_energy - law.rhs) == 0)

def print():
    return pretty(law, use_unicode=False)

@validate_input(energy_=units.energy, time_=units.time)
@validate_output(units.power)
def calculate_power(energy_: Quantity, time_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({energy: energy_, time: time_})
    return expr_to_quantity(result_expr, "power")