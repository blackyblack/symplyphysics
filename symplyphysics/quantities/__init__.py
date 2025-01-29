"""
Physical constants
==================

Contains useful physical constants. Fundamental constants can be also found in `sympy.physics.units`_ module.

.. _sympy.physics.units: https://github.com/sympy/sympy/blob/b6376e085a047b8bada988ff57fbb79a62546842/sympy/physics/units/definitions/unit_definitions.py#L245
"""

from sympy.physics import units
from symplyphysics.core.symbols.quantities import Quantity

standard_conditions_temperature = Quantity(273.15 * units.kelvin,
    display_symbol="t_std",
    display_latex="t_\\text{std}")
"""
Zero Celsius degrees. The :symbols:`temperature` at which water freezes.
It is also temperature for Standard Temperature and Pressure (STP)

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Standard_temperature_and_pressure>`__.
"""

standard_laboratory_temperature = Quantity(298 * units.kelvin,
    display_symbol="t_lab",
    display_latex="t_\\text{lab}")
"""
Approximately :math:`25` degrees Celsius. Commonly used :symbols:`temperature` for tabulation purposes.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Room_temperature#Definitions_in_science_and_industry>`__.
"""

electron_rest_mass = Quantity(9.1093837015e-31 * units.kilogram,
    display_symbol="m_e",
    display_latex="m_\\text{e}")
"""
:symbols:`mass` of stationary electron.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electron_mass>`__.
"""

bohr_radius = Quantity(0.529e-10 * units.meter, display_symbol="a0", display_latex="a_0")
"""
The Bohr radius is the radius of the electron orbit of the hydrogen atom closest
to the nucleus in the atomic model proposed by Niels Bohr.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Bohr_radius>`__.
"""

hydrogen_ionization_energy = Quantity(13.6 * units.electronvolt,
    display_symbol="IE_h",
    display_latex="\\mathrm{IE}_\\text{H}")
"""
The ionization energy is the smallest energy required to remove an electron from
a free atom in its basic energy state to infinity.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Ionization_energy#Bohr_model_for_hydrogen_atom>`__.
"""

solar_mass = Quantity(1.9884e30 * units.kilogram, display_symbol="M_Sun", display_latex="M_\\odot")
r"""
The solar :symbols:`mass` is a standard unit of mass in astronomy approximately equal to the mass
of the Sun. The relative uncertainty of the measurement is :math:`4 \cdot 10^{-5}`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Solar_mass>`__.
"""

earth_mass = Quantity(5.9722e24 * units.kilogram,
    display_symbol="M_Earth",
    display_latex="M_\\oplus")
"""
The Earth :symbols:`mass` is a standard unit of mass in astronomy equal to the mass of the planet 
Earth. The relative uncertainty of the measurement is :math:`10^{-4}`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Earth_mass>`__.
"""

boltzmann_constant = Quantity(units.boltzmann_constant,
    display_symbol="k_B",
    display_latex="k_\\text{B}")
"""
The Boltzmann constant is the proportionality factor that relates the average relative thermal energy of particles
in a gas with the thermodynamic temperature of the gas.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Boltzmann_constant>`__.
"""

molar_gas_constant = Quantity(units.molar_gas_constant, display_symbol="R")
"""
The gas constant is the constant of proportionality that relates the energy scale in physics to the temperature
scale and the scale used for amount of substance.
It is molar equivalent to the :attr:`~boltzmann_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Gas_constant>`__.
"""

speed_of_light = Quantity(units.speed_of_light, display_symbol="c")
r"""
The speed of light in vacuum is a universal physical constant that is exactly equal to :math:`299 \, 792 \, 458` metres per second.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Speed_of_light>`__.
"""

vacuum_permittivity = Quantity(
    units.vacuum_permittivity,
    display_symbol="epsilon_0",
    display_latex="\\varepsilon_0",
)
"""
Vacuum permittivity, also known as permittivity of free space or the electric constant, is the value of the absolute
dielectric permittivity of classical vacuum.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Vacuum_permittivity>`__.
"""

vacuum_permeability = Quantity(units.vacuum_permeability,
    display_symbol="mu_0",
    display_latex="\\mu_0")
"""
Vacuum permeability, also known as permeability of free space of the magnetic constant, is the value of the absolute
permeability of classical vacuum.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Vacuum_permeability>`__.
"""

elementary_charge = Quantity(units.elementary_charge, display_symbol="e")
"""
**Elementary charge** is a fundamental physical constant defined as the electric charge carried by a single proton
or, equivalently, the magnitude of the negative charge carried by a single electron.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Elementary_charge>`__.
"""

hbar = Quantity(units.hbar, display_symbol="hbar", display_latex="\\hbar")
"""
**Reduced Planck constant** is a modified version of the Planck constant used in the description of
Quantum Mechanics.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Planck_constant#Reduced_Planck_constant>`__.
"""

planck = Quantity(units.planck, display_symbol="h")
"""
The **Planck constant** is a fundamental physical constant of foundational importance in quantum mechanics.
It is the constant of proportionality between a photon's energy and frequency.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Planck_constant#>`__.
"""

avogadro_constant = Quantity(units.avogadro, display_symbol="N_A", display_latex="N_\\text{A}")
"""
The **Avogadro constant** is an SI defining constant defined as the number of constituent particles
per mole.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Avogadro_constant>`__.
"""

acceleration_due_to_gravity = Quantity(units.acceleration_due_to_gravity, display_symbol="g")
"""
A conventional standard value of the gravitational acceleration at Earth's surface.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Standard_gravity>`__.
"""

stefan_boltzmann_constant = Quantity(units.stefan_boltzmann_constant,
    display_symbol="sigma",
    display_latex="\\sigma")
"""
The **Stefan—Boltzmann constant** is the constant of proportionality between radiant exitance and black body's temperature
in the Stefan—Boltzmann law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Stefan%E2%80%93Boltzmann_law#>`__.
"""

richardson_constant = Quantity(120.17 * (units.ampere / units.kelvin**2 / units.centimeter**2),
    display_symbol="a")
"""
Constant of proportionality proposed by Richardson to describe the law of thermionic emission.

**Links:**

#. `Richardson's law <https://en.wikipedia.org/wiki/Thermionic_emission#Richardson's_law>`__.
"""

rydberg_frequency = Quantity(3.2898419602500e15 * units.hertz,
    display_symbol="R_H",
    display_latex="R_\\text{H}")
"""
In spectroscopy, the **Rydberg constant** is a physical constant relating to the electromagnetic spectra of an atom.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Rydberg_constant#Rydberg_frequency>`__.
"""

wien_displacement_constant = Quantity(
    (speed_of_light * planck) / (boltzmann_constant * 4.965114),
    display_symbol="b",
)
"""
A constant of proportionality in *Wien's displacement law*.

**Links:**

#. `Wien's displacement law <https://en.wikipedia.org/wiki/Wien%27s_displacement_law>`__.
"""

gravitational_constant = Quantity(units.gravitational_constant, display_symbol="G")
"""
The **gravitational constant** is a physical constant used in calculating the gravitational
attraction between two objects.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Gravitational_constant>`__.
"""

hubble_constant = Quantity(7e-11 / units.year, display_symbol="H")
"""
The **Hubble's constant** is the proportionality constant between the recessional velocity and
the proper distance between the galaxy and the observer in the Hubble's law. Its exact value
is up to debate, however, the fundamental theory gives the number :math:`7%` per billion years.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Hubble%27s_law>`__.
"""

zero_point_luminosity = Quantity(3.0128e28 * units.watt, display_symbol="L_0")
"""
**Zero-point luminosity** is a constant defined relative to a star for calibrating perposes. The
value given here has been defined by the International Astronomical Union (IAU).

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Zero_point_(photometry)#Bolometric_magnitude_zero_point>`__.
"""

sun_luminosity = Quantity(3.827e26 * units.watt, display_symbol="L_Sun", display_latex="L_\\odot")
"""
:symbols:`luminosity` of the Sun.
"""

faraday_constant = Quantity(elementary_charge * avogadro_constant,
    display_symbol="F",
    display_latex="\\mathfrak{F}")
"""
The **Faraday constant** represents the amount of electric charge carried by one mole,
or Avogadro's number, of electrons.

**Links:**

#. `TechTarget <https://www.techtarget.com/whatis/definition/Faraday-constant>`__.
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
    "rydberg_frequency",
    "wien_displacement_constant",
    "hubble_constant",
    "zero_point_luminosity",
    "faraday_constant",
]
