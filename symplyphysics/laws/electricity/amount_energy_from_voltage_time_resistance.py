from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity,
)
# Description

# The amount of energy released by a conductor with a current is directly proportional
# to the square of the applied voltage, the time of the current and inversely proportional
# to the resistance of the conductor. This is the differential form of the law
## This Joule-Lenz law is derived from more generic law:
## See: [dissipated_heat_power](./dissipated_heat_power_is_proportional_to_current_square.py) implementation.

# Amount of energy Q = U**2 * t / R
# where:
# U - voltage to conductor
# t - time of current action
# R - resistance of conductor

amount_energy, voltage, time, resistance = symbols(
    'amount_energy voltage time resistance')
law = Eq(amount_energy, (voltage**2 * time) / resistance)

def print():
    return pretty(law, use_unicode=False)

@validate_input(voltage_=units.voltage, time_=units.time, resistance_=units.impedance)
@validate_output(units.energy)
def calculate_amount_energy(voltage_: Quantity, time_: Quantity, resistance_: Quantity) -> Quantity:
    result_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_energy_expr.subs({
        voltage: voltage_, time: time_, resistance: resistance_})
    return expr_to_quantity(result_expr, 'amount_energy')