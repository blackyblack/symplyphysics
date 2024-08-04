#!/usr/bin/env python3

from sympy import symbols, Eq, solve, S
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, units
from symplyphysics.laws.thermodynamics.fermi_dirac_statistics import (
    single_particle_state_distribution as distribution_law,)

occupancy = distribution_law.occupancy_of_state
energy = distribution_law.energy_of_state
chemical_potential = distribution_law.total_chemical_potential
temperature = distribution_law.temperature

# Expressing energy and temperature in the units of chemical potential

reduced_energy = symbols("reduced_energy")
reduced_energy_eqn = Eq(reduced_energy, energy / chemical_potential)

reduced_temperature = symbols("reduced_temperature")
reduced_temperature_eqn = Eq(
    reduced_temperature,
    units.boltzmann_constant * temperature / chemical_potential,
)

distribution_expr = solve(
    [distribution_law.law, reduced_energy_eqn, reduced_temperature_eqn],
    (occupancy, energy, temperature),
    dict=True,
)[0][occupancy]

print(f"Occupancy as a function of reduced energy:\n\n{print_expression(distribution_expr)}\n")

base_plot = plot(
    title="Occupancy as a function of reduced energy for Fermi—Dirac distribution",
    xlabel=r"$\varepsilon/\mu$",
    ylabel=r"$\bar n$",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)

factors_ = 1, 2, 10, 100

for factor_ in factors_:
    expr_ = distribution_expr.subs(reduced_temperature, S.One / factor_)
    label = r"$k_\text{B} T = \mu" + (f"/{factor_}$" if factor_ != 1 else "$")
    sub_plot = plot(
        expr_,
        (reduced_energy, 0, 5),
        label=label,
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
