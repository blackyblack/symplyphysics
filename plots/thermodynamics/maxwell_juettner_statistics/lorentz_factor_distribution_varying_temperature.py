#!/usr/bin/env python3

from sympy import solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, units
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.thermodynamics.relativistic import (
    maxwell_juettner_distribution as distribution_law,
    reduced_temperature_in_maxwell_juettner_statistics as reduced_law,
)

print("Maxwell-Juettner distribution formula:\n")
print(print_expression(distribution_law.law))

law = distribution_law.law.subs(
    distribution_law.reduced_temperature,
    reduced_law.reduced_temperature,
)

rhs = solve(
    (law, reduced_law.law),
    (distribution_law.distribution_function, reduced_law.reduced_temperature),
    dict=True,
)[0][distribution_law.distribution_function]

rhs = rhs.subs(reduced_law.particle_mass, units.electron_rest_mass)

rhs = evaluate_expression(rhs)

temperatures_ = [1e10, 2e10, 5e10]  # K

base_plot = plot(
    title="Maxwell—Jüttner distribution of Lorentz factor at different temperatures",
    xlabel=r"Lorentz factor $\gamma$",
    ylabel=r"Probability $f(\gamma)$",
    legend=True,
    backend=MatplotlibBackend,
    show=False,
)

for temperature_ in temperatures_:
    expr = rhs.subs(
        reduced_law.temperature,
        temperature_,
    )
    sub_plot = plot(
        expr,
        (distribution_law.lorentz_factor, 1, 100),
        label=f"T = {temperature_:.2e} K",
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
