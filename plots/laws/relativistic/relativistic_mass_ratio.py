#!/usr/bin/env python3

from sympy import solve, symbols, Eq
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from sympy.physics.units import speed_of_light
import symplyphysics.laws.relativistic.relativistic_mass as relativistic_mass

# Description
## Plot the dependency of the relativistic to rest mass ratio on the fraction of the speed of light

mass_ratio, speed_ratio = symbols("mass_ratio beta", positive=True)
rel_mass, rest_mass = symbols("rel_mass rest_mass", positive=True)
speed = symbols("speed", positive=True)

# plain definition of ratio
mass_ratio_definition = Eq(mass_ratio, rel_mass / rest_mass)
speed_ratio_definition = Eq(speed_ratio, speed / speed_of_light)

rel_mass_sub = solve(mass_ratio_definition, rel_mass)[0]
speed_sub = solve(speed_ratio_definition, speed)[0]

law = relativistic_mass.law.subs({
    relativistic_mass.relativistic_mass: rel_mass_sub,
    relativistic_mass.rest_mass: rest_mass,
    relativistic_mass.velocity: speed_sub,
})

result_mass_ratio = solve(law, mass_ratio)[0]
solved_law = Eq(mass_ratio, result_mass_ratio)

print(f"Formula is:\n{print_expression(solved_law)}")

p1 = plot(
    result_mass_ratio,
    (speed_ratio, 0, 1),
    ylim=(0, 100),
    axis_center=(0.0, 0.0),
    line_color="red",
    title="Relativistic to rest mass ratio depending on speed ratio",
    xlabel=r"$\frac{v}{c}$",
    ylabel=r"$\frac{m_{rel}}{m_{rest}}$",
    legend=True,
    backend=MatplotlibBackend,
    show=False,
)
p1.show()
