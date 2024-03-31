from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

# Description
## Potential energy of body EP = m * g * h
## Where:
## m - body mass
## h - height from Earth surface
## g - free fall acceleration

potential_energy_of_body = Symbol("potential_energy_of_body", units.energy)
height = Symbol("height", units.length)
free_fall_acceleration = units.acceleration_due_to_gravity

law = Eq(potential_energy_of_body, symbols.basic.mass * free_fall_acceleration * height)


def print_law() -> str:
    return print_expression(law)


@validate_input(body_mass_=symbols.basic.mass, height_=height)
@validate_output(potential_energy_of_body)
def calculate_potential_energy(body_mass_: Quantity, height_: Quantity) -> Quantity:
    result_energy_expr = solve(law, potential_energy_of_body,
        dict=True)[0][potential_energy_of_body]
    result_expr = result_energy_expr.subs({symbols.basic.mass: body_mass_, height: height_})
    return Quantity(result_expr)
