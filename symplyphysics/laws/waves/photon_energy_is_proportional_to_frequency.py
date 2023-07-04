from sympy import (Eq, solve)
from sympy.physics.units import planck as planck_constant
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output)

# Description
## Photon is the elementary part of any electromagnetical radiation which has no mass and always moves with speed of light.
## It carries some energy. The amount of energy depends only on the frequency of the photon.
## Law is: E = h * nu, where
## E is energy of photon,
## h is Planck constant,
## nu is frequency of photon.

photon_energy = Symbol("photon_energy", units.energy)
photon_frequency = Symbol("frequency", units.frequency)

law = Eq(photon_energy, planck_constant * photon_frequency)


def print_law() -> str:
    return print_expression(law)


@validate_input(photon_frequency_=photon_frequency)
@validate_output(photon_energy)
def calculate_energy(photon_frequency_: Quantity) -> Quantity:
    result_energy_expr = solve(law, photon_energy, dict=True)[0][photon_energy]
    result_expr = result_energy_expr.subs({photon_frequency: photon_frequency_})
    return expr_to_quantity(result_expr)
