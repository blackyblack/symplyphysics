#!/usr/bin/env python3

from sympy import Idx, solve, Symbol, Eq
from symplyphysics import print_expression, Quantity, units, convert_to, global_index
from symplyphysics.definitions import mass_fraction as mass_fraction_law
from symplyphysics.laws.conservation import mass_after_equals_to_mass_before as conservation_mass_law
from symplyphysics.laws.conservation import mixture_mass_equal_sum_of_components_masses as mixture_mass_law

# Example from 6th in https://uchitel.pro/%D0%B7%D0%B0%D0%B4%D0%B0%D1%87%D0%B8-%D0%BD%D0%B0-%D0%B2%D1%8B%D1%87%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%B8%D0%B5-%D0%BC%D0%B0%D1%81%D1%81%D1%8B-%D1%80%D0%B0%D1%81%D1%82%D0%B2%D0%BE%D1%80%D0%B5%D0%BD/
# The mass fraction of salt in a solution saturated at a temperature of 40 ° C is 35%.
# When 300 g of this solution is cooled to a temperature of 10 ° C, 45 g of salt precipitates.
# Determine the solubility of the salt at 10 ° C.

mass_of_sediment = Symbol("mass_of_sediment")
mass_fraction_of_salt = Symbol("mass_fraction_of_salt")
mass_of_mixture = Symbol("mass_of_mixture")

mass_of_water = Symbol("mass_of_water")
mass_of_salt_after = Symbol("mass_of_salt_after")

solubility = Symbol("solubility")

mass_fraction_equation = mass_fraction_law.definition.subs({
    mass_fraction_law.mass_of_mixture: mass_of_mixture,
    mass_fraction_law.mass_fraction: mass_fraction_of_salt
})
mass_of_salt_before_value = solve(mass_fraction_equation,
    mass_fraction_law.mass_of_component,
    dict=True)[0][mass_fraction_law.mass_of_component]

index_local = Idx("index_local", (1, 2))
mass_of_two_components = mixture_mass_law.law.subs(global_index, index_local).doit()

mass_of_salt_and_sediment = mass_of_two_components.subs({
    mixture_mass_law.component_mass[1]: mass_of_salt_after,
    mixture_mass_law.component_mass[2]: mass_of_sediment,
}).rhs
conservation_salt_mass_equation = conservation_mass_law.law.subs({
    conservation_mass_law.mass(conservation_mass_law.initial_time):
        mass_of_salt_before_value,
    conservation_mass_law.mass(conservation_mass_law.final_time):
        mass_of_salt_and_sediment
})
mass_of_salt_after_value = solve(conservation_salt_mass_equation, mass_of_salt_after,
    dict=True)[0][mass_of_salt_after]

mass_of_mixture_equation = mass_of_two_components.subs({
    mixture_mass_law.component_mass[1]: mass_of_water,
    mixture_mass_law.component_mass[2]: mass_of_salt_before_value,
    mixture_mass_law.mixture_mass: mass_of_mixture,
})

mass_of_water_value = solve(mass_of_mixture_equation, mass_of_water, dict=True)[0][mass_of_water]

solubility_value = mass_fraction_law.definition.subs({
    mass_fraction_law.mass_of_mixture: mass_of_water_value,
    mass_fraction_law.mass_of_component: mass_of_salt_after_value
}).rhs
answer = Eq(solubility, solubility_value)
print(f"Total equation:\n{print_expression(answer)}")

solubility_fraction = solubility_value.subs({
    mass_of_mixture: Quantity(300 * units.grams),
    mass_of_sediment: Quantity(45 * units.grams),
    mass_fraction_of_salt: Quantity(35 * units.percents)
})
print(
    f"Solubility is: {print_expression(convert_to(Quantity(solubility_fraction), units.percents).evalf(4))} %"
)
