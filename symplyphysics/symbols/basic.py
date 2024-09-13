"""
Basic
=====

Symbols of fundamental physical quantities.
"""

from symplyphysics import units, SymbolNew

time = SymbolNew("t", units.time)
"""
**Time** is a scalar physical quantity operationally defined as the reading of a clock, specifically
a count of repeating events such as the SI second. It is a fundamental concept used to define other
physical quantities.
"""

period = SymbolNew("T", units.time)
"""
**Period** is the duration of time of one cycle in a repeating event.
"""

mass = SymbolNew("m", units.mass)
"""
**Mass** is an intrinsic scalar property of a body, and one can distinguish at least seven different aspects
of mass that define it. Experiments have shown that these values are proportional, and in some cases
equal. Often, the inertial mass is being used, which is a measure of the object's resistance to acceleration
when a force is applied.
"""

work = SymbolNew("W", units.energy)
"""
**Work** is the energy transferred to or from an object via the application of force.
"""

energy_density = SymbolNew("w", units.energy / units.volume)
"""
Energy per unit volume.
"""

energy = SymbolNew("E", units.energy)
"""
**Energy** is the quantitative property that is transferred to a body or to a physical system, recognizable in
the performance of work and in the form of heat and light.
"""

power = SymbolNew("P", units.power)
"""
**Power** is the amount of energy transferred or converted per unit time.
"""

radius_of_curvature = SymbolNew("r", units.length)
"""
**Radius of curvature** is the inverse of curvature and is equal to the distance to the center of curvature.
"""

density = SymbolNew("rho", units.mass / units.volume, display_latex="\\rho")
"""
**Density** is mass per unit volume.
"""
