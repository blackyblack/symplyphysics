from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol,
                           print_expression, validate_input_symbols,
                           validate_output_symbol)

# Description
# Kinetic energy of body: EK = (m * v**2) / 2
# Where:
# m - body mass
# v - body velocity

kinetic_energy_of_body = Symbol("kinetic_energy_of_body", units.energy)
body_mass = Symbol("body_mass", units.mass)
body_velocity = Symbol("body_velocity", units.velocity)

law = Eq(kinetic_energy_of_body, body_mass * body_velocity**2 / 2)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(body_mass_=body_mass, body_velocity_=body_velocity)
@validate_output_symbol(kinetic_energy_of_body)
def calculate_kinetic_energy(body_mass_: Quantity,
                             body_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy_of_body,
                               dict=True)[0][kinetic_energy_of_body]
    result_expr = result_energy_expr.subs({
        body_mass: body_mass_,
        body_velocity: body_velocity_
    })
    return expr_to_quantity(result_expr)
