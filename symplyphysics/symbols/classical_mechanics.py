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

position = SymbolNew("x", units.length)
"""
**Position** is defined as the spatial location of an object with respect to a coordinate system.
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

thickness = SymbolNew("d", units.length)
"""
**Thickness** is a measure of a size of an object, usually the separation between two layers, or
the distance through an object distinct from length and width.
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

phase_speed = SymbolNew("v", units.speed)
"""
**Phase speed** is the speed at which the phase of the wave travels.
"""

pressure = SymbolNew("p", units.pressure)
"""
**Pressure** is the force applied perpendicular to the surface of an object per unit area over which
that force is distributed.
"""

temporal_frequency = SymbolNew("f", units.frequency)
"""
**Temporal frequency** is the number of occurrences of a repeating event per unit of time.
"""

sound_intensity_level = SymbolNew("L_I", dimensionless)
"""
**Sound intensity level** is the measure of the *intensity* of a sound relative to a reference value.
"""

rotational_inertia = SymbolNew("I", units.mass * units.length**2)
"""
**Rotational inertia**, also known as **moment of inertia**, is defined relative to a rotational axis
and is the ratio between the torque applied and the resulting angular acceleration about that axis.
"""

quality_factor = SymbolNew("Q", dimensionless)
"""
**Quality factor** or **Q factor** is a dimensionless parameter that describes how underdamped an
oscillator or resonator is.
"""

momentum = SymbolNew("p", units.momentum)
"""
**Momentum**, more specifically **linear** or **transitional momentum**, is the product of the mass and
velocity of an object.
"""

mechanical_energy = SymbolNew("E", units.energy)
"""
**Mechanical energy** is defined to be the sum of potential energy and kinetic energy.
"""

kinetic_energy = SymbolNew("K", units.energy)
"""
**Kinetic energy** of an object is the form of energy that it possesses due to its motion.
"""

potential_energy = SymbolNew("U", units.energy)
"""
**Potential energy** is the energy held by an object because of its position relative to other objects,
stresses within itself, its electric charge, or other factors. Potential energy is associated with so
called conservative forces and only depends on the initial and final positions of the body in space.
"""

mass_flow_rate = SymbolNew("mu", units.mass / units.time, display_latex="\\mu")
"""
**Mass flow rate** is the mass of a substance which passes per unit of time.
"""

stiffness = SymbolNew("k", units.force / units.length)
"""
**Stiffness** is the extent to which an object resists deformation in response to an applied force.
"""

compliance = SymbolNew("c", units.length / units.force)
"""
**Compliance** is the inverse of :symbols:`stiffness`.
"""
