"""
Thermodynamics (Symbols)
========================

Symbols related to thermodynamics.
"""

from sympy.physics import units
from ..core.symbols.symbols import SymbolNew
from ..core.dimensions import dimensionless

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

thermal_expansion_coefficient = SymbolNew("alpha", 1 / units.temperature, display_latex="\\alpha")
"""
**Thermal expansion coefficient** describes how the size of an object changes with a change in temperature at
constant pressure.
"""

thermodynamic_compressibility = SymbolNew("beta", 1 / units.pressure, display_latex="\\beta")
"""
**Compressibility** is a measure of the instantaneous relative volume change of a fluid or solid as a response
to a pressure or mean stress change.
"""
