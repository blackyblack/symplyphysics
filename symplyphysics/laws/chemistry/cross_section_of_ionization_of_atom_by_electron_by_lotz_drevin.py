from sympy import Eq, solve, pi, log
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, dimensionless)
from symplyphysics.quantities import bohr_radius, hydrogen_ionization_energy

# Description
## The effective cross section is a physical quantity characterizing the probability of transition of a system of
## two interacting particles to a certain final state, a quantitative characteristic of the acts of collision of
## particles of a stream hitting a target with target particles. The effective cross-section has the dimension of the area.
## In this law, we are talking about the interaction of an atom and an electron, which ionizes an atom.
## Equivalent electrons on the outer shell of the ionized atom are electrons with the same principal and orbital quantum numbers.
## In this case, the Lotz-Drewin approximation for the ionization cross section is considered.

## Law is: g = 2.66 * pi * a0^2 * l * (Uh^2 / Ui^2) * B1 * ((U / Ui - 1) / (U / Ui)^2) * ln(1.25 * B2 * U / Ui), where
## g - the cross-sectional area of the ionization of particles,
## Ui - the ionization energy of atom,
## U - the energy of ionizing electrons,
## B1 - first calculation coefficient,
## B2- second calculation coefficient,
## l - the number of equivalent electrons on the outer shell of the ionized atom,
## Uh - hydrogen ionization energy,
## a0 - Bohr radius.

cross_sectional_area_of_ionization = Symbol("cross_sectional_area_of_ionization", units.area)

ionization_energy = Symbol("ionization_energy", units.energy)
energy_of_electron = Symbol("energy_of_electron", units.energy)
first_calculation_coefficient = Symbol("first_calculation_coefficient", dimensionless)
second_calculation_coefficient = Symbol("second_calculation_coefficient", dimensionless)
number_of_equivalent_electrons_on_outer_orbit = Symbol(
    "number_of_equivalent_electrons_on_outer_orbit", dimensionless)

expression_1 = energy_of_electron / ionization_energy
expression_2 = 2.66 * pi * bohr_radius**2 * number_of_equivalent_electrons_on_outer_orbit * hydrogen_ionization_energy**2 / ionization_energy**2
expression_3 = first_calculation_coefficient * (expression_1 - 1) / expression_1**2
expression_4 = log(1.25 * second_calculation_coefficient * expression_1)
law = Eq(cross_sectional_area_of_ionization, expression_2 * expression_3 * expression_4)


@validate_input(ionization_energy_=ionization_energy,
    energy_of_electron_=energy_of_electron,
    first_calculation_coefficient_=first_calculation_coefficient,
    second_calculation_coefficient_=second_calculation_coefficient,
    number_of_equivalent_electrons_on_outer_orbit_=number_of_equivalent_electrons_on_outer_orbit)
@validate_output(cross_sectional_area_of_ionization)
def calculate_cross_sectional_area_of_ionization(
        ionization_energy_: Quantity, energy_of_electron_: Quantity,
        first_calculation_coefficient_: float, second_calculation_coefficient_: float,
        number_of_equivalent_electrons_on_outer_orbit_: int) -> Quantity:
    result_expr = solve(law, cross_sectional_area_of_ionization,
        dict=True)[0][cross_sectional_area_of_ionization]
    result_expr = result_expr.subs({
        ionization_energy:
            ionization_energy_,
        energy_of_electron:
            energy_of_electron_,
        first_calculation_coefficient:
            first_calculation_coefficient_,
        second_calculation_coefficient:
            second_calculation_coefficient_,
        number_of_equivalent_electrons_on_outer_orbit:
            number_of_equivalent_electrons_on_outer_orbit_,
    })
    return Quantity(result_expr)
