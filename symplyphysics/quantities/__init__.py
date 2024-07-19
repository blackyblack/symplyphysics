"""
Physical constants
==================

Contains useful physical constants. Fundamental constants can be also found in `sympy.physics.units`_ module.

.. _sympy.physics.units: https://github.com/sympy/sympy/blob/b6376e085a047b8bada988ff57fbb79a62546842/sympy/physics/units/definitions/unit_definitions.py#L245
"""

from sympy.physics import units
from ..core.symbols.quantities import Quantity

standard_conditions_temperature = Quantity(273.15 * units.kelvin)
r"""
Zero Celsius degrees. The temperature at which water freezes.
It is also temperature for Standard Temperature and Pressure (STP)

Symbol:
    t_std

Latex:
    :math:`t_\text{std}`
"""

standard_laboratory_temperature = Quantity(298 * units.kelvin)
r"""
Approximately 25 Celsius degrees. Commonly used temperature for tabulation purposes.

Symbol:
    t_lab

Latex:
    :math:`t_\text{lab}`
"""

electron_rest_mass = Quantity(9.1093837015e-31 * units.kilogram)
r"""
Mass of stationary electron

Symbol:
    m_e

Latex:
    :math:`m_\text{e}`
"""

bohr_radius = Quantity(0.529e-10 * units.meter)
r"""
The Bohr radius is the radius of the electron orbit of the hydrogen atom closest
to the nucleus in the atomic model proposed by Niels Bohr.

Symbol:
    a0

Latex:
    :math:`a_0`
"""

hydrogen_ionization_energy = Quantity(13.6 * units.electronvolt)
r"""
The ionization energy is the smallest energy required to remove an electron from
a free atom in its basic energy state to infinity.

Symbol:
    IE_h

Latex:
    :math:`\mathrm{IE}_\text{H}`
"""

solar_mass = Quantity(1.9884e30 * units.kilogram)
r"""
The solar mass is a standard unit of mass in astronomy approximately equal to the mass
of the Sun.

Symbol:
    M_Sun

Latex:
    :math:`\mathrm{M}_\odot`
"""

__all__ = [
    "standard_conditions_temperature",
    "standard_laboratory_temperature",
    "electron_rest_mass",
    "bohr_radius",
    "hydrogen_ionization_energy",
    "solar_mass",
]
