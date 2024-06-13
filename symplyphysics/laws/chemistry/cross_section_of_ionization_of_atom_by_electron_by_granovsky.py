from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The effective cross section is a physical quantity characterizing the probability of transition of a system of
## two interacting particles to a certain final state, a quantitative characteristic of the acts of collision of
## particles of a stream hitting a target with target particles. The effective cross-section has the dimension of the area.
## In this law, we are talking about the interaction of an atom and an electron, which ionizes an atom.
## In this case, the Granovsky approximation for the ionization cross section is considered.
## As the electron energy increases, the velocities of primary and secondary electrons increase, the possibility of their
## recombination with ions decreases, and the ionization cross-section area increases. However, at very high energies, the
## ionization cross section decreases, as the electrons fly past the atom without having time to ionize it, since the time spent
## by the electron near the atom decreases. Therefore, there is a maximum ionization cross section.

## Law is: g = g_max * ((U - Ui) / (Umax - Ui)) * exp((Umax - U) / (Umax - Ui)), where
## g - the cross-sectional area of the ionization of particles,
## g_max - the maximum value of the ionization cross-sectional area,
## Ui - the ionization energy of atom,
## U - the energy of ionizing electrons,
## Umax - the electron energy corresponding to the maximum value of the ionization cross-sectional area.

cross_sectional_area_of_ionization = Symbol("cross_sectional_area_of_ionization", units.area)

ionization_energy = Symbol("ionization_energy", units.energy)
energy_of_electron = Symbol("energy_of_electron", units.energy)
maximum_cross_sectional_area_of_ionization = Symbol("maximum_cross_sectional_area_of_ionization", units.area)
energy_of_electron_at_max_area = Symbol("energy_of_electron_at_max_area", units.energy)

expression_1 = (energy_of_electron - ionization_energy) / (energy_of_electron_at_max_area - ionization_energy)
expression_2 = maximum_cross_sectional_area_of_ionization * expression_1
expression_3 = exp((energy_of_electron_at_max_area - energy_of_electron) / (energy_of_electron_at_max_area - ionization_energy))

law = Eq(cross_sectional_area_of_ionization, expression_2 * expression_3)


@validate_input(ionization_energy_=ionization_energy,
    energy_of_electron_=energy_of_electron,
    maximum_cross_sectional_area_of_ionization_=maximum_cross_sectional_area_of_ionization,
    energy_of_electron_at_max_area_=energy_of_electron_at_max_area)
@validate_output(cross_sectional_area_of_ionization)
def calculate_cross_sectional_area_of_ionization(ionization_energy_: Quantity,
    energy_of_electron_: Quantity, maximum_cross_sectional_area_of_ionization_: Quantity,
    energy_of_electron_at_max_area_: Quantity) -> Quantity:
    result_expr = solve(law, cross_sectional_area_of_ionization,
        dict=True)[0][cross_sectional_area_of_ionization]
    result_expr = result_expr.subs({
        ionization_energy: ionization_energy_,
        energy_of_electron: energy_of_electron_,
        maximum_cross_sectional_area_of_ionization: maximum_cross_sectional_area_of_ionization_,
        energy_of_electron_at_max_area: energy_of_electron_at_max_area_,
    })
    return Quantity(result_expr)
