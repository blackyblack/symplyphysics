from sympy import Eq, solve
from sympy.physics.units import speed_of_light as c
from symplyphysics import units, Quantity, Symbol, print_expression, validate_input, validate_output

# Description
## Fundamentally inner energy of an object is synonimical to its mass.

# Law: E = m * c**2, where
## E is rest energy of body/system,
## m is mass,
## c is speed of light.

rest_energy = Symbol("rest_energy", units.energy)
rest_mass = Symbol("rest_mass", units.mass)

law = Eq(rest_energy, rest_mass * c**2)


def print_law() -> str:
    return print_expression(law)


@validate_input(rest_mass_=rest_mass)
@validate_output(rest_energy)
def calculate_rest_energy(rest_mass_: Quantity) -> Quantity:
    result_expr = solve(law, rest_energy, dict=True)[0][rest_energy]
    energy_applied = result_expr.subs({rest_mass: rest_mass_})
    return Quantity(energy_applied)
