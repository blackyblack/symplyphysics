#!/usr/bin/env python3

from sympy import symbols, Eq, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, units
from symplyphysics.laws.thermodynamics.bose_einstein_statistics import (
    single_particle_state_distribution as distribution_law,
)

occupancy = distribution_law.occupancy_of_state
energy = distribution_law.energy_of_state
chemical_potential = distribution_law.total_chemical_potential
temperature = distribution_law.temperature

reduced_energy = symbols("reduced_energy")
reduced_energy_eqn = Eq(reduced_energy, energy / (units.boltzmann_constant * temperature))

reduced_chemical_potential = symbols("reduced_chemical_potential")
reduced_chemical_potential_eqn = Eq(
    reduced_chemical_potential,
    chemical_potential / (units.boltzmann_constant * temperature)
)

distribution_expr = solve(
    [distribution_law.law, reduced_energy_eqn, reduced_chemical_potential_eqn],
    (occupancy, energy, chemical_potential),
    dict=True,
)[0][occupancy]

print(f"Occupancy as a function of reduced energy and chemical potential:\n\n{print_expression(distribution_expr)}\n")

base_plot = plot(
    title="Occupancy as a function of reduced energy and chemical potential for Bose—Eintein distribution",
    xlabel=r"$-\varepsilon / k_\text{B} T$",
    ylabel=r"$\bar n$",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)

values_ = -0.5, -1, -2, -5

for value_ in values_:
    expr_ = distribution_expr.subs(reduced_chemical_potential, value_)
    label = r"$\mu = " + (str(value_) if value_ != -1 else "-") + r" k_\text{B} T$"
    sub_plot = plot(
        expr_,
        (reduced_energy, 0, 5),
        label=label,
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
