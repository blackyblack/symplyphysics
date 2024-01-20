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
## Electric potential (also called the electric field potential, potential drop, the electrostatic potential)
## is defined as the amount of work energy needed per unit of electric charge to move the charge from
## a reference point to a specific point in an electric field. More precisely, the electric potential is
## the energy per unit charge for a test charge that is so small that the disturbance of the field under
## consideration is negligible.

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
