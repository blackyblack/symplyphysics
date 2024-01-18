from sympy import (Eq, solve, pi)
from sympy.physics.units import electric_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    dimensionless,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Potential energy is a function of the state of the system. The zero value is taken when the charges are infinitely
## far from each other. Note also that this is the energy of the entire system, the energy of interaction, so it makes
## no sense to say that some part of this energy belongs to one of the charges.

## Law is: W = q1 * q2 / (4 * pi  * e0 * e * r), where
## W - energy of interaction of charges,
## q1 - first charge,
## q2 - second charge,
## e0 - electric constant,
## e - relative permittivity of medium,
## r - distance between charges.

energy = Symbol("energy", units.energy)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
distance = Symbol("distance", units.length)
charge_1 = Symbol("charge_1", units.charge)
charge_2 = Symbol("charge_2", units.charge)

law = Eq(
    energy,
    charge_1 * charge_2 / (4 * pi * electric_constant * relative_permittivity * distance))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity, distance_=distance, charge_1_=charge_1, charge_2_=charge_2)
@validate_output(energy)
def calculate_energy(relative_permittivity_: float, distance_: Quantity, charge_1_: Quantity, charge_2_: Quantity) -> Quantity:
    result_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        distance: distance_,
        charge_1: charge_1_,
        charge_2: charge_2_,
    })
    return Quantity(result_expr)
