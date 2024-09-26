#!/usr/bin/env python3
from sympy import solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.thermodynamics import radiance_of_black_body_from_temperature as stefan_boltzmann_law

print(f"Formula is:\n{print_expression(stefan_boltzmann_law.law)}")

solved = solve(stefan_boltzmann_law.law, stefan_boltzmann_law.radiance,
    dict=True)[0][stefan_boltzmann_law.radiance]

radiance_temperature = evaluate_expression(solved)
print(f"Radiance function is:\n{print_expression(radiance_temperature)}")

p1 = plot(radiance_temperature, (stefan_boltzmann_law.temperature, 0, 400),
    line_color="black",
    title="Stefan - Boltzmann Law",
    xlabel="T(K)",
    ylabel="P(W/m^2)",
    label="completely black body",
    backend=MatplotlibBackend,
    legend=True,
    show=False)
p1.show()
