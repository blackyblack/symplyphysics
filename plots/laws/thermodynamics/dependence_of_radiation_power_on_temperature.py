#!/usr/bin/env python3
from sympy import solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.laws.thermodynamics import stefan_boltzmann_law

print(f"Formula is:\n{stefan_boltzmann_law.print_law()}")

temperature_function = stefan_boltzmann_law.temperature

solved = solve(stefan_boltzmann_law.law, stefan_boltzmann_law.irradiance,
    dict=True)[0][stefan_boltzmann_law.irradiance]
irradiance_temperature = solved.subs({
    stefan_boltzmann_law.temperature: temperature_function,
    stefan_boltzmann_law.units.stefan: 5.67036713*(10**(-8)),
})

print(f"Pressure function is:\n{print_expression(irradiance_temperature)}")

p1 = plot(irradiance_temperature, (temperature_function, 0, 400),
    line_color="black",
    title="Stefan - Boltzmann Law",
    xlabel="T(K)",
    ylabel="P(W/m^2)",
    label="completely black body",
    backend=MatplotlibBackend,
    legend=True,
    show=False)
p1.show()
