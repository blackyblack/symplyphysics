#!/usr/bin/env python3

from sympy import symbols, Eq, solve, S
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, units
from symplyphysics.laws.thermodynamics.fermi_dirac_statistics import (
    single_particle_state_distribution as distribution_law,
)

occupancy = distribution_law.occupancy_of_state
energy = distribution_law.energy_of_state
chemical_potential = distribution_law.total_chemical_potential
temperature = distribution_law.temperature

reduced_temperature = symbols("reduced_temperature")
reduced_temperature_eqn = Eq(
    reduced_temperature,
    units.boltzmann_constant * temperature / (energy - chemical_potential),
)

distribution_expr = solve(
    [distribution_law.law, reduced_temperature_eqn],
    (occupancy, temperature),
    dict=True
)[0][occupancy]

print(f"Occupancy as a function of reduced temperature:\n{print_expression(distribution_expr)}\n")

base_plot = plot(
    distribution_expr,
    (reduced_temperature, 0, 15),
    title="Occupancy as a function of reduced temperature for Fermi—Dirac distribution",
    xlabel=r"$\frac{k_\text{B} T}{\varepsilon - \mu}$",
    ylabel=r"$\bar n$",
    backend=MatplotlibBackend,
    show=False,
)

asymptote_expr = distribution_expr.limit(reduced_temperature, S.Infinity)

asymptote_plot = plot(
    asymptote_expr,
    (reduced_temperature, 0, 15),
    show=False,
)
base_plot.extend(asymptote_plot)

base_plot.show()
