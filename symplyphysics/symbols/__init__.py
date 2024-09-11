"""
Physical symbols
================

Symbols represent physical quantities, units, mathematical operations and relationships.
"""

from symplyphysics import units, SymbolNew

mass = SymbolNew("m", units.mass)
"""
Mass is an intrinsic scalar property of a body, and one can distinguish at least seven different aspects
of mass that define it. Experiments have shown that these values are proportional, and in some cases
equal. Often, the inertial mass is being used, which is a measure of the object's resistance to acceleration
when a force is applied.
"""

force = SymbolNew("F", units.force)
"""
Force is a vector quantity denoting an influence that can cause an object to change its velocity
unless counterbalanced by other forces.
"""

acceleration = SymbolNew("a", units.acceleration)
"""
Acceleration is the rate of change of the object's velocity with respect to time. It is a vector quantity.
"""

temperature = SymbolNew("T", units.temperature)
"""
Temperature is a scalar quantity that quantitatively expresses the attribute of hotness and coldness. It reflects
the average kinetic energy of the vibrating and colliding atoms making up a substance.
"""

__all__ = [
    "mass",
    "force",
    "acceleration",
    "temperature",
]
