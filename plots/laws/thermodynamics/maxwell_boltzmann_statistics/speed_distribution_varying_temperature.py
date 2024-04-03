#!/usr/bin/env python3

from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, convert_to, Quantity, units
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

argon_mass_ = convert_to(Quantity(39.948 * units.amu), units.kilogram)
temperatures_ = [100, 200, 300, 400, 500]  # K

print(f"Maxwell-Boltzmann speed distribution function of Argon:\n{print_expression(speed_distribution.law)}\n")

temperature_plot = plot(
    title="Maxwell—Boltzmann speed distribution of Argon at different temperatures",
    xlabel=r"speed $v, \frac{m}{s}$",
    ylabel=r"probability density, $(\frac{m}{s})^{-1}$",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)

distribution = speed_distribution.law.rhs.subs({
    speed_distribution.particle_mass: argon_mass_,
    units.boltzmann_constant: convert_to(units.boltzmann_constant, units.joule / units.kelvin),
})

for temperature_ in temperatures_:
    temperature_subplot = plot(
        distribution.subs(speed_distribution.equilibrium_temperature, temperature_),
        (speed_distribution.particle_speed, 0, 1000),
        label=f"$T = {temperature_} K$",
        show=False,
    )
    temperature_plot.append(temperature_subplot[0])

temperature_plot.show()
