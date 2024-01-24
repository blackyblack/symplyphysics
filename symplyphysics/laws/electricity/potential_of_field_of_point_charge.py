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
## The potential of the electrostatic field φ at a given point is a scalar value equal to the ratio
## of the potential energy W of a charge q placed at a given point to the value of this charge.
## The unit of measurement of potential is – Volt (V). A volt is equal to the potential of the field point
## at which a charge of 1 Coulomb has a potential energy of 1 Joule.

## Law is: φ = q / (4 * pi * e0 * e * r), where
## φ - potential of electrostatic field,
## q - charge,
## e0 - electric constant,
## e - electric relative permittivity,
## r - distance from point charge.

electrostatic_potential = Symbol("electrostatic_potential", units.voltage)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
distance = Symbol("distance", units.length)
charge = Symbol("charge", units.charge)

law = Eq(electrostatic_potential,
    charge / (4 * pi * electric_constant * relative_permittivity * distance))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity, distance_=distance, charge_=charge)
@validate_output(electrostatic_potential)
def calculate_electrostatic_potential(relative_permittivity_: float, distance_: Quantity,
    charge_: Quantity) -> Quantity:
    result_expr = solve(law, electrostatic_potential, dict=True)[0][electrostatic_potential]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        distance: distance_,
        charge: charge_,
    })
    return Quantity(result_expr)
