from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Capacitor can accumlate energy in the electric field inside it.
## Law: Q = C * U^2 / 2
## Q is accumulated energy
## C is capacitance of capacitor
## U is voltage on the capacitor

accumulated_energy, capacitance, voltage = symbols('accumulated_energy capacitance voltage')
law = Eq(accumulated_energy, capacitance * voltage**2 / 2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(capacitance_=units.capacitance, voltage_=units.voltage)
@validate_output(units.energy)
def calculate_accumulated_energy(capacitance_: Quantity, voltage_: Quantity) -> Quantity:
    result_energy_expr = solve(law, accumulated_energy, dict=True)[0][accumulated_energy]
    result_expr = result_energy_expr.subs({capacitance: capacitance_, voltage: voltage_})
    return expr_to_quantity(result_expr, 'accumulated_energy')
