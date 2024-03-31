#!/usr/bin/env python3
from sympy import solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (print_expression, units, Quantity, convert_to)
from symplyphysics.laws.thermodynamics import radiance_of_black_body_from_temperature as stefan_boltzmann_law

print(f"Formula is:\n{stefan_boltzmann_law.print_law()}")

solved = solve(stefan_boltzmann_law.law, stefan_boltzmann_law.radiance,
    dict=True)[0][stefan_boltzmann_law.radiance]

radiance_temperature = solved.subs({
    stefan_boltzmann_law.units.stefan_boltzmann_constant:
    convert_to(Quantity(units.stefan_boltzmann_constant),
    units.watt / units.meter**2 / units.kelvin**4).evalf(6)
})
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
