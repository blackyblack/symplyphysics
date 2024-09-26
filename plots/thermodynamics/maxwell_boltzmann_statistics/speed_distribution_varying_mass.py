#!/usr/bin/env python3

from collections import namedtuple
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, Quantity, units, quantities, convert_to_si
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

MassDatum = namedtuple("MassDatum", "mass label")

# Used [this python library](https://pymatgen.org/pymatgen.core.html#pymatgen.core.composition.Composition)
# to get the mass data in atomic mass units (amu).
mass_data_ = [
    MassDatum(mass=4.0026, label=r"\text{He}"),
    MassDatum(mass=17.031, label=r"\text{NH}_3"),
    MassDatum(mass=28.010, label=r"\text{CO}"),
    MassDatum(mass=39.948, label=r"\text{Ar}"),
    MassDatum(mass=44.096, label=r"\text{C}_3\text{H}_8"),
]

print(
    f"Maxwell-Boltzmann speed distribution function:\n{print_expression(speed_distribution.law)}\n")

mass_plot = plot(
    title="Maxwell—Boltzmann speed distribution for particles of different masses, $T$ = const",
    xlabel=r"speed $v, \frac{m}{s}$",
    ylabel=r"probability density, $(\frac{m}{s})^{-1}$",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)

distribution = speed_distribution.law.rhs.subs(
    speed_distribution.equilibrium_temperature, quantities.standard_conditions_temperature
)
distribution = evaluate_expression(distribution)

for mass_datum_ in mass_data_:
    mass_ = convert_to_si(Quantity(mass_datum_.mass * units.amu))
    mass_subplot = plot(
        distribution.subs(speed_distribution.particle_mass, mass_),
        (speed_distribution.particle_speed, 0, 2000),
        label=f"$m_{{{mass_datum_.label}}} = {mass_datum_.mass} \\, \\text{{amu}}$",
        show=False,
    )
    mass_plot.append(mass_subplot[0])

mass_plot.show()
