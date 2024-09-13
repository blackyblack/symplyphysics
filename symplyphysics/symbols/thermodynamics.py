"""
Thermodynamics
==============

Symbols related to thermodynamics.
"""

from symplyphysics import units, SymbolNew, dimensionless

temperature = SymbolNew("T", units.temperature)
"""
Temperature is a scalar quantity that quantitatively expresses the attribute of hotness and coldness. It reflects
the average kinetic energy of the vibrating and colliding atoms making up a substance.
"""

adiabatic_index = SymbolNew("gamma", dimensionless, display_latex="\\gamma")
"""
**Adiabatic index**, or **heat capacity ratio**, is the ratio of heat capacity at constant pressure to that
at constant volume.
"""

heat_capacity = SymbolNew("C", units.energy / units.temperature)
"""
**Heat capacity** or **thermal capacity** is a physical property of matter, defined as the amount of heat to be
supplied to an object to produce a unit change in its temperature.
"""
