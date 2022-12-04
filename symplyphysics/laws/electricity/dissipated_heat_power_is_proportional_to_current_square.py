from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, simplify, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.electricity import power_from_energy_time as power_and_time
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law
from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law

# Description
# Dissipated heat power is proportional to current square and resistance
# P = I**2 * R
# where :
# P - is dissipated heat power of the element
# I is current flowing through this element
# R is impedance of this element

heat_power, current, resistance = symbols('heat_power current resistance')
law = Eq(heat_power, current**2 * resistance)

# This law might be easily derived via Joule-Lenz law and dependence power from energy and time

voltage_applied = solve(ohm_law.law, ohm_law.voltage, dict=True)[0][ohm_law.voltage]
energy_applied =  solve(power_and_time.law, power_and_time.energy, dict=True)[0][power_and_time.energy]
law_applied = joule_lenz_law.law.subs({joule_lenz_law.amount_energy: energy_applied,
    joule_lenz_law.voltage: voltage_applied
})
law_derived = solve(law_applied, power_and_time.power, dict=True)[0][power_and_time.power]

# Check if derived power is same as declared
difference = simplify(law_derived - law.rhs)
assert(difference == 0)

def print():
    return pretty(law_derived, use_unicode=False)

@validate_input(current_=units.current, resistance_=units.impedance)
@validate_output(units.power)
def calculate_heat_power(current_: Quantity, resistance_: Quantity) -> Quantity:
    result_power_expr = solve(law, heat_power, dict=True)[0][heat_power]
    result_expr = result_power_expr.subs({current: current_, resistance: resistance_})
    return expr_to_quantity(result_expr, 'heat_power')
