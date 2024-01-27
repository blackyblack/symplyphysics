from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The electrostatic potential is a scalar physical quantity equal to the ratio of the potential energy
## of the interaction of a charge with an electric field to the magnitude of the charge.
## The potential energy is equal to the work that a fixed charge would do on a test charge if
## the test charge moved from the initial distance "r" to infinity.

## Law is: ф = W / q, where
## ф - electrostatic potential,
## W - potential energy,
## q - charge.

electrostatic_potential = Symbol("electrostatic_potential", units.voltage)

potential_energy = Symbol("potential_energy", units.energy)
charge = Symbol("charge", units.charge)

law = Eq(electrostatic_potential, potential_energy / charge)


def print_law() -> str:
    return print_expression(law)


@validate_input(potential_energy_=potential_energy, charge_=charge)
@validate_output(electrostatic_potential)
def calculate_potential(potential_energy_: Quantity, charge_: Quantity) -> Quantity:
    result_expr = solve(law, electrostatic_potential, dict=True)[0][electrostatic_potential]
    result_expr = result_expr.subs({
        potential_energy: potential_energy_,
        charge: charge_,
    })
    return Quantity(result_expr)
