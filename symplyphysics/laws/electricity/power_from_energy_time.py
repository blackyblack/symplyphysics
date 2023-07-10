from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
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

power = Symbol("power", units.power)
energy = Symbol("energy", units.energy)
time = Symbol("time", units.time)

law = Eq(power, energy / time)

# Derive the same law from definition of power as energy derivative

linear_function_coefficient = Symbol("linear_function_coefficient", units.energy / units.time)
initial_energy_constant = Symbol("initial_energy_constant", units.energy)

energy_linear_function = linear_function_coefficient * time + initial_energy_constant

power_definition_applied = power_derivative.definition.subs(power_derivative.time, time)
power_definition_applied = power_definition_applied.subs({
    power_derivative.energy(time): energy_linear_function,
    power_derivative.power(time): power
})
power_applied_eq = power_definition_applied.doit()
energy_eq = Eq(energy, energy_linear_function)

# derived_power = (energy - initial_energy_constant)/time
derived_power = solve([power_applied_eq, energy_eq], (power, linear_function_coefficient),
    dict=True)[0][power]

# Assume initial_energy_constant = 0, according to law conditions
derived_power_without_initial_energy = derived_power.subs(initial_energy_constant, 0)

# Check that derived power is same as declared
assert expr_equals(derived_power_without_initial_energy, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(energy_=energy, time_=time)
@validate_output(power)
def calculate_power(energy_: Quantity, time_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({energy: energy_, time: time_})
    return Quantity(result_expr)
