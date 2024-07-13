"""
Kinematics
==========

Contains properties for 'kinematics' branch of physics

Kinematics is a subfield of physics and mathematics, developed in classical mechanics, that describes the motion of points, bodies (objects), and systems of bodies (groups of objects) without considering the forces that cause them to move.
"""

from sympy.physics import units
from ...core.symbols.symbols import Symbol

acceleration = Symbol("acceleration", units.acceleration)
"""
Acceleration is the rate of change of the velocity of an object with respect to time.

Symbols:
    a
"""

__all__ = [
    "acceleration",
]
