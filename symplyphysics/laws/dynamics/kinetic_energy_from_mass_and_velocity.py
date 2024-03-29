from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

# Description
# Kinetic energy of body: EK = (m * v**2) / 2
# Where:
# m - body mass
# v - body velocity

kinetic_energy_of_body = Symbol("kinetic_energy_of_body", units.energy)
body_velocity = Symbol("body_velocity", units.velocity)

law = Eq(kinetic_energy_of_body, symbols.basic.mass * body_velocity**2 / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(body_mass_=symbols.basic.mass, body_velocity_=body_velocity)
@validate_output(kinetic_energy_of_body)
def calculate_kinetic_energy(body_mass_: Quantity, body_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy_of_body, dict=True)[0][kinetic_energy_of_body]
    result_expr = result_energy_expr.subs({
        symbols.basic.mass: body_mass_,
        body_velocity: body_velocity_
    })
    return Quantity(result_expr)
