"""
Classical mechanics
===================

Symbols related to classical mechanics.
"""

from symplyphysics import units, SymbolNew

force = SymbolNew("F", units.force)
"""
Force is a vector quantity denoting an influence that can cause an object to change its velocity
unless counterbalanced by other forces.
"""

speed = SymbolNew("v", units.velocity)
"""
Speed is the rate of change of the object's position with respect to time. It is the scalar counterpart
of the velocity vector, which is the rate of change of the object's displacement with respect to time.
"""

acceleration = SymbolNew("a", units.acceleration)
"""
Acceleration is the rate of change of the object's velocity with respect to time. It is a vector quantity.
"""

distance = SymbolNew("d", units.length)
"""
Distance is is a measure of the spatial separation between two points.
"""

radial_distance = SymbolNew("r", units.length)
"""
Distance to the origin of the coordinate system.
"""
