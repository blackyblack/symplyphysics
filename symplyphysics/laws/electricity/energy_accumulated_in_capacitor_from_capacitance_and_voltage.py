from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol,
                           print_expression, validate_input_symbols,
                           validate_output_symbol)

# Description
## Capacitor can accumlate energy in the electric field inside it.
## Law: Q = C * U^2 / 2
## Q is accumulated energy
## C is capacitance of capacitor
## U is voltage on the capacitor

accumulated_energy = Symbol("accumulated_energy", units.energy)
capacitance = Symbol("capacitance", units.capacitance)
voltage = Symbol("voltage", units.voltage)

law = Eq(accumulated_energy, capacitance * voltage**2 / 2)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(capacitance_=capacitance, voltage_=voltage)
@validate_output_symbol(accumulated_energy)
def calculate_accumulated_energy(capacitance_: Quantity,
                                 voltage_: Quantity) -> Quantity:
    result_energy_expr = solve(law, accumulated_energy,
                               dict=True)[0][accumulated_energy]
    result_expr = result_energy_expr.subs({
        capacitance: capacitance_,
        voltage: voltage_
    })
    return expr_to_quantity(result_expr)
