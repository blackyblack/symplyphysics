from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.definitions import power_from_energy_time as power_and_time
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law
from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law

# Description
## Dissipated heat power is proportional to current square and resistance
## Joule-Lenz law: P = I**2 * R
## P - is dissipated heat power of the element
## I is current flowing through this element
## R is impedance of this element

power = joule_lenz_law.law.subs({joule_lenz_law.amount_energy: power_and_time.energy, joule_lenz_law.voltage: ohm_law.voltage })
law = solve(power, power_and_time.energy, dict=True)[0][power_and_time.energy]

def print():
    return pretty(law, use_unicode=False)

@validate_input(current_=units.current, resistance_=units.impedance)
@validate_output(units.power)
def calculate_heat_power(current_: Quantity, resistance_: Quantity) -> Quantity:
    result_power_expr = solve(law, heat_power, dict=True)[0][heat_power]
    result_expr = result_power_expr.subs({current: current_, resistance: resistance_})
    return expr_to_quantity(result_expr, 'heat_power')
