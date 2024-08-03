#!/usr/bin/env python3
from sympy import sqrt, symbols, Eq, solve, limit, oo
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, units
from symplyphysics.laws.condensed_matter import (
    concentration_of_intrinsic_charge_carriers as carriers_law,)

reduced_concentration = symbols("reduced_concentration")

density_product_expr = (carriers_law.density_of_states_in_conduction_band *
    carriers_law.density_of_states_in_valence_band)

reduced_concentration_eqn = Eq(
    reduced_concentration,
    carriers_law.charge_carriers_concentration / sqrt(density_product_expr),
)

reduced_temperature = symbols("reduced_temperature")

reduced_temperature_eqn = Eq(
    reduced_temperature,
    (units.boltzmann_constant * carriers_law.temperature) / carriers_law.band_gap,
)

reduced_concentration_expr = solve(
    (carriers_law.law, reduced_concentration_eqn, reduced_temperature_eqn),
    (carriers_law.charge_carriers_concentration, reduced_concentration, carriers_law.temperature),
    dict=True,
)[0][reduced_concentration]

reduced_concentration_limit = limit(reduced_concentration_expr, reduced_temperature, oo)

print("Definition of reduced concentration:")
print(print_expression(reduced_concentration_eqn))

print("\nDefinition of reduced temperature gap:")
print(print_expression(reduced_temperature_eqn))

print("\nThe formula for reduced concentration of charge carriers:\n")
print(print_expression(reduced_concentration_expr))

base_plot = plot(
    title="Reduced concentration of intrinsic charge carriers versus reduced temperature",
    xlabel=r"reduced temperature $T^* = \frac{k_\text{B} T}{E_g}$",
    ylabel=r"reduced concentration $n^* = \frac{n}{\sqrt{N_c N_v}}$",
    legend=True,
    backend=MatplotlibBackend,
    show=False,
)

main_plot = plot(
    reduced_concentration_expr,
    (reduced_temperature, 0, 5),
    label=r"$n^*(T^*)$",
    show=False,
)
base_plot.extend(main_plot)

limit_plot = plot(
    reduced_concentration_limit,
    (reduced_temperature, 0, 5),
    label=r"$n^*$ as $T^* \to \infty$",
    show=False,
)
base_plot.extend(limit_plot)

base_plot.show()
