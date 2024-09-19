"""
Classical mechanics (Symbols)
=============================

Symbols related to classical mechanics.
"""

from sympy.physics import units
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from ..core.symbols.symbols import SymbolNew
from ..core.dimensions import dimensionless

force = SymbolNew("F", units.force)
"""
**Force** is a vector quantity denoting an influence that can cause an object to change its velocity
unless counterbalanced by other forces.
"""

speed = SymbolNew("v", units.velocity)
"""
**Speed** is the rate of change of the object's position with respect to time. It is the scalar counterpart
of the velocity vector, which is the rate of change of the object's displacement with respect to time.
"""

acceleration = SymbolNew("a", units.acceleration)
"""
**Acceleration** is the rate of change of the object's velocity with respect to time. It is a vector quantity.
"""

distance = SymbolNew("d", units.length)
"""
**Distance** is is a measure of the spatial separation between two points.
"""

distance_to_origin = SymbolNew("r", units.length)
"""
Distance to the origin of the coordinate system.
"""

distance_to_axis = SymbolNew("r", units.length)
"""
Distance to reference axis.
"""

length = SymbolNew("l", units.length)
"""
**Length** is a measure of a size of an object.
"""

area = SymbolNew("A", units.area)
"""
**Area** is the size of a region on a two-dimensional surface.
"""

angular_speed = SymbolNew("w", angle_type / units.time, display_latex="\\omega")
"""
**Angular speed** is the rate of change of angular distance with respect to time.
"""

angular_frequency = angular_speed
"""
**Angular frequency**, also called **angular rate**, is a scalar measure of the temporal
rate of change of the phase argument of a sinusoidal waveform or sine function.
"""

angular_acceleration = SymbolNew("alpha", angle_type / units.time**2, display_latex="\\alpha")
"""
**Angular acceleration** is the rate of change of angular speed with respect to time.
"""

angular_distance = SymbolNew("theta", angle_type, display_latex="\\theta")
"""
**Angular distance** is a measure of an angular separation between two points.
"""

angular_wavenumber = SymbolNew("k", angle_type / units.length)
"""
**Angular wavenumber** is the spatial analog of temporal frequency equal to radians per unit length.
"""

wavelength = SymbolNew("lambda", units.length, display_latex="\\lambda")
"""
**Wavelength** or **spatial period** is the distance over which the wave's shape repeats.
"""

damping_ratio = SymbolNew("zeta", dimensionless, display_latex="\\zeta")
"""
**Damping ratio** is a dimensionless measure describing how oscillations in a system decay after a disturbance.
"""

volume = SymbolNew("V", units.volume)
"""
Volume is a measure of regions in three-dimensional space.
"""

impulse = SymbolNew("J", units.momentum)
"""
**Impulse** is the change in momentum of an object.
"""
