"""
Physical symbols
================

Symbols represent physical quantities, units, mathematical operations and relationships.
"""

from symplyphysics import units, SymbolNew

# General

time = SymbolNew("t", units.time)
"""
Time is a scalar physical quantity operationally defined as the reading of a clock, specifically
a count of repeating events such as the SI second. It is a fundamental concept used to define other
physical quantities.
"""

mass = SymbolNew("m", units.mass)
"""
Mass is an intrinsic scalar property of a body, and one can distinguish at least seven different aspects
of mass that define it. Experiments have shown that these values are proportional, and in some cases
equal. Often, the inertial mass is being used, which is a measure of the object's resistance to acceleration
when a force is applied.
"""

# Dynamics

force = SymbolNew("F", units.force)
"""
Force is a vector quantity denoting an influence that can cause an object to change its velocity
unless counterbalanced by other forces.
"""

# Kinematics

speed = SymbolNew("v", units.velocity)
"""
Speed is the rate of change of the object's position with respect to time. It is the scalar counterpart
of the velocity vector, which is the rate of change of the object's displacement with respect to time.
"""

acceleration = SymbolNew("a", units.acceleration)
"""
Acceleration is the rate of change of the object's velocity with respect to time. It is a vector quantity.
"""

# Thermodynamics

temperature = SymbolNew("T", units.temperature)
"""
Temperature is a scalar quantity that quantitatively expresses the attribute of hotness and coldness. It reflects
the average kinetic energy of the vibrating and colliding atoms making up a substance.
"""

# Electrodynamics

admittance = SymbolNew("Y", units.conductance)
"""
Admittance is a measure of how easily a circuit or device will allow a current to flow, defined as the reciprocal
of impedance.
"""

electrical_impedance = SymbolNew("Z", units.impedance)
"""
Electrical impedance is the opposition to current presented by the combined effect of
resistance and reactance in a circuit.
"""


__all__ = [
    "time",
    "mass",
    "force",
    "speed",
    "acceleration",
    "temperature",
]
