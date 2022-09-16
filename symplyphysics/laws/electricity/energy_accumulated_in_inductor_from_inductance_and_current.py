from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Inductor can accumlate energy in the magnetic field inside it.
## Law: Q = L * I^2 / 2
## Q is accumulated energy
## L is inductance of inductor
## I is current flowing through the inductor

accumulated_energy, inductance, current = symbols('accumulated_energy inductance current')
law = Eq(accumulated_energy, inductance * current**2 / 2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(inductance_=units.inductance, current_=units.current)
@validate_output(units.energy)
def calculate_accumulated_energy(inductance_: Quantity, current_: Quantity) -> Quantity:
    result_energy_expr = solve(law, accumulated_energy, dict=True)[0][accumulated_energy]
    result_expr = result_energy_expr.subs({inductance: inductance_, current: current_})
    return expr_to_quantity(result_expr, 'accumulated_energy')
