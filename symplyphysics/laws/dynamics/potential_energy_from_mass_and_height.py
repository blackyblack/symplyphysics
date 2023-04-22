from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol,
                           print_expression, validate_input_symbols,
                           validate_output_symbol)

# Description
## Potential energy of body EP = m * g * h
## Where:
## m - body mass
## h - height from Earth surface
## g - free fall acceleration

potential_energy_of_body = Symbol("potential_energy_of_body", units.energy)
height = Symbol("height", units.length)
body_mass = Symbol("body_mass", units.mass)
free_fall_acceleration = units.acceleration_due_to_gravity

law = Eq(potential_energy_of_body, body_mass * free_fall_acceleration * height)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(body_mass_=body_mass, height_=height)
@validate_output_symbol(potential_energy_of_body)
def calculate_potential_energy(body_mass_: Quantity,
                               height_: Quantity) -> Quantity:
    result_energy_expr = solve(law, potential_energy_of_body,
                               dict=True)[0][potential_energy_of_body]
    result_expr = result_energy_expr.subs({
        body_mass: body_mass_,
        height: height_
    })
    return expr_to_quantity(result_expr)
