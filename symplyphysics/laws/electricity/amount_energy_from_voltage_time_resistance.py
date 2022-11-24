from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity,
)
# Description
# Amount heat energy of current Q = U^2 *t / R
# where:
# U - voltage
# t - time
# R - resistance
amount_energy, operating_voltage, operating_time, resistance_heater = symbols(
    'amount_energy, operating_voltage, operating_time, resistance_heater')
law = Eq(amount_energy, (operating_voltage**2 * operating_time) / resistance_heater)

def print():
    return pretty(law, use_unicode=False)

@validate_input(operating_voltage_=units.voltage, operating_time_=units.time, resistance_heater_=units.impedance)
@validate_output(units.energy)
def calculate_amount_energy(operating_voltage_: Quantity, operating_time_: Quantity, resistance_heater_: Quantity) -> Quantity:
    result_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_energy_expr.subs({
        operating_voltage: operating_voltage_, operating_time: operating_time_, resistance_heater: resistance_heater_})
    return expr_to_quantity(result_expr, 'amount_energy')