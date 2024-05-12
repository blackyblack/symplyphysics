from sympy import Eq, solve, pi
from sympy.physics.units import electric_constant, elementary_charge
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output,)

# Description
## The effective cross section is a physical quantity characterizing the probability of transition of a system of
## two interacting particles to a certain final state, a quantitative characteristic of the acts of collision of
## particles of a stream hitting a target with target particles. The effective cross-section has the dimension of the area.
## In a magnetron, this value can be calculated if you know the ionization energy of gas atoms.

## Law is: g =  q^2 / (2 * pi * e0^2 * Ui^2), where
## g - the cross-sectional area of the interaction of particles,
## q - elementary charge,
## e0 - electric constant,
## Ui - the ionization energy of atoms, expressed in volts.

cross_sectional_area_of_interaction = Symbol("cross_sectional_area_of_interaction", units.area)

ionization_energy = Symbol("ionization_energy", units.voltage)

law = Eq(cross_sectional_area_of_interaction, elementary_charge**2 / (2 * pi * electric_constant**2 * ionization_energy**2))


@validate_input(ionization_energy_=ionization_energy)
@validate_output(cross_sectional_area_of_interaction)
def calculate_cross_sectional_area_of_interaction(ionization_energy_: Quantity) -> Quantity:
    result_expr = solve(law, cross_sectional_area_of_interaction, dict=True)[0][cross_sectional_area_of_interaction]
    result_expr = result_expr.subs({
        ionization_energy: ionization_energy_,
    })
    return Quantity(result_expr)
