#!/usr/bin/env python3
from sympy import solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (print_expression, units, convert_to)
from symplyphysics.laws.thermodynamics import irradiation_of_black_body_from_temperature as stefan_boltzmann_law

print(f"Formula is:\n{stefan_boltzmann_law.print_law()}")
solved = solve(stefan_boltzmann_law.law, stefan_boltzmann_law.irradiance,
    dict=True)[0][stefan_boltzmann_law.irradiance]
irradiance_temperature = solved.subs({
    stefan_boltzmann_law.units.stefan_boltzmann_constant: 
                                convert_to(units.stefan_boltzmann_constant,units.watt/units.meter**2/units.kelvin**4).evalf(5)
})

print(f"Pressure function is:\n{print_expression(irradiance_temperature_)}")
p1 = plot(irradiance_temperature, (stefan_boltzmann_law.temperature, 0, 400),   
    line_color="black",
    title="Stefan - Boltzmann Law",
    xlabel="T(K)",
    ylabel="P(W/m^2)",
    label="completely black body",
    backend=MatplotlibBackend,
    legend=True,
    show=False)
p1.show()
