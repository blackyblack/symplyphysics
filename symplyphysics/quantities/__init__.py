"""
Physical constants
==================

Contains useful physical constants. Fundamental constants can be also found in `sympy.physics.units`_ module.

.. _sympy.physics.units: https://github.com/sympy/sympy/blob/b6376e085a047b8bada988ff57fbb79a62546842/sympy/physics/units/definitions/unit_definitions.py#L245
"""

from sympy.physics import units
from ..core.symbols.quantities import Quantity

standard_conditions_temperature = Quantity(273.15 * units.kelvin,
    display_symbol="t_std",
    display_latex="t_\\text{std}")
"""
Zero Celsius degrees. The :attr:`~symplyphysics.symbols.temperature` at which water freezes.
It is also temperature for Standard Temperature and Pressure (STP)
"""

standard_laboratory_temperature = Quantity(298 * units.kelvin,
    display_symbol="t_lab",
    display_latex="t_\\text{lab}")
"""
Approximately 25 Celsius degrees. Commonly used :attr:`~symplyphysics.symbols.temperature` for tabulation purposes.
"""

electron_rest_mass = Quantity(9.1093837015e-31 * units.kilogram,
    display_symbol="m_e",
    display_latex="m_\\text{e}")
"""
:attr:`~symplyphysics.symbols.mass` of stationary electron
"""

bohr_radius = Quantity(0.529e-10 * units.meter, display_symbol="a0", display_latex="a_0")
"""
The Bohr radius is the radius of the electron orbit of the hydrogen atom closest
to the nucleus in the atomic model proposed by Niels Bohr.
"""

hydrogen_ionization_energy = Quantity(13.6 * units.electronvolt,
    display_symbol="IE_h",
    display_latex="\\mathrm{IE}_\\text{H}")
"""
The ionization energy is the smallest energy required to remove an electron from
a free atom in its basic energy state to infinity.
"""

solar_mass = Quantity(1.9884e30 * units.kilogram, display_symbol="M_Sun", display_latex="M_\\odot")
r"""
The solar :attr:`~symplyphysics.symbols.mass` is a standard unit of mass in astronomy approximately equal to the mass
of the Sun. The relative uncertainty of the measurement is :math:`4 \cdot 10^{-5}`.
"""

earth_mass = Quantity(5.9722e24 * units.kilogram,
    display_symbol="M_Earth",
    display_latex="M_\\oplus")
"""
The Earth :attr:`~symplyphysics.symbols.mass` is a standard unit of mass in astronomy equal to the mass of the planet 
Earth. The relative uncertainty of the measurement is :math:`10^{-4}`.
"""

boltzmann_constant = Quantity(units.boltzmann_constant, display_symbol="k_B", display_latex="k_\\text{B}")
"""
The Boltzmann constant is the proportionality factor that relates the average relative thermal energy of particles
in a gas with the thermodynamic temperature of the gas.
"""

molar_gas_constant = Quantity(units.molar_gas_constant, display_symbol="R")
"""
The gas constant is the constant of proportionality that relates the energy scale in physics to the temperature
scale and the scale used for amount of substance.
It is molar equivalent to the :attr:`~symplyphysics.quantities.boltzmann_constant`
"""

speed_of_light = Quantity(units.speed_of_light, display_symbol="c")
"""
The speed of light in vacuum is a universal physical constant that is exactly equal to 299,792,458 metres per second.
"""

vacuum_permittivity = Quantity(units.vacuum_permittivity,
    display_symbol="epsilon_0",
    display_latex="\\varepsilon_0",
)
"""
Vacuum permittivity, also known as permittivity of free space or the electric constant, is the value of the absolute
dielectric permittivity of classical vacuum.
"""

vacuum_permeability = Quantity(units.vacuum_permeability,
    display_symbol="mu_0",
    display_latex="\\mu_0")
"""
Vacuum permeability, also known as permeability of free space of the magnetic constant, is the value of the absolute
permeability of classical vacuum.
"""

elementary_charge = Quantity(units.elementary_charge, display_symbol="e")
"""
**Elementary charge** is a fundamental physical constant defined as the electric charge carried by a single proton
or, equivalently, the magnitude of the negative charge carried by a single electron.
"""

hbar = Quantity(units.hbar, display_symbol="hbar", display_latex="\\hbar")
"""
**Reduced Planck constant** is a modified version of the Planck constant used in the description of
Quantum Mechanics.
"""

planck = Quantity(units.planck, display_symbol="h")
"""
The **Planck constant** is a fundamental physical constant of foundational importance in quantum mechanics.
It is the constant of proportionality between a photon's energy and frequency.
"""

avogadro_constant = Quantity(units.avogadro, display_symbol="N_A", display_latex="N_\\text{A}")
"""
The **Avogadro constant** is an SI defining constant defined as the number of constituent particles
per mole.
"""

acceleration_due_to_gravity = Quantity(units.acceleration_due_to_gravity, display_symbol="g")
"""
A conventional standard value of the gravitational acceleration at Earth's surface.
"""

stefan_boltzmann_constant = Quantity(units.stefan_boltzmann_constant, display_symbol="sigma", display_latex="\\sigma")
"""
The **Stefan—Boltzmann constant** is the constant of proportionality between radiant exitance and black body's temperature
in the Stefan—Boltzmann law.
"""

richardson_constant = Quantity(
    120.17 * (units.ampere / units.kelvin**2 / units.centimeter**2),
    display_symbol="a")
"""
Constant of proportionality proposed by Richardson to describe the law of thermionic emission.

**Links:**

#. `Richardson's law <https://en.wikipedia.org/wiki/Thermionic_emission#Richardson's_law>`__.
"""


__all__ = [
    "standard_conditions_temperature",
    "standard_laboratory_temperature",
    "electron_rest_mass",
    "bohr_radius",
    "hydrogen_ionization_energy",
    "solar_mass",
    "earth_mass",
    "boltzmann_constant",
    "molar_gas_constant",
    "speed_of_light",
    "vacuum_permittivity",
    "vacuum_permeability",
    "elementary_charge",
    "hbar",
    "planck",
    "avogadro_constant",
    "acceleration_due_to_gravity",
    "stefan_boltzmann_constant",
    "richardson_constant",
]
