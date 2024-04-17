#!/usr/bin/env python3

from sympy import symbols, Eq, solve, S
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

# Expressing temperature in reduced units using the quantities in the law
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
    title="Occupancy as a function of reduced temperature for Bose-Einstein distribution",
    xlabel=r"reduced temperature $T^* = \frac{k_\text{B} T}{\varepsilon - \mu}$",
    ylabel=r"occupancy $\bar n$",
    legend=True,
    backend=MatplotlibBackend,
    show=False,
)

occupancy_plot = plot(
    distribution_expr,
    (reduced_temperature, 0, 15),
    label=r"$\bar n (T^*)$",
    show=False,
)
base_plot.extend(occupancy_plot)

asymptote_expr = distribution_expr.series(reduced_temperature, S.Infinity, 1).removeO()

asymptote_plot = plot(
    asymptote_expr,
    (reduced_temperature, 0, 15),
    label=r"$T^* \to \infty$ asymptote",
    show=False,
)
base_plot.extend(asymptote_plot)

base_plot.show()

