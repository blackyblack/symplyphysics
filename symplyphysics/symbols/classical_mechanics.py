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

young_modulus = SymbolNew("E", units.pressure)
"""
**Young modulus** is a mechanical property of solid materials that measures the tensile or compressive
stiffness when the force is applied lengthwise.
"""

phase_shift = SymbolNew("phi", angle_type, display_latex="\\varphi")
"""
**Phase shift**, also known as **phase offset** or **phase difference**, is the shift of phase between
two periodic functions.
"""

coefficient_of_friction = SymbolNew("mu", dimensionless, display_latex="\\mu")
"""
**Coefficient of friction** is a dimensionless scalar value which equals the ratio of the force of
friction between two bodies and the force pressing them together, either during or at the onset of
slipping.
"""

height = SymbolNew("h", units.length)
"""
**Height** is measure of vertical distance, either vertical extent or vertical position.

In the case of three-dimensional space, height is measured along the vertical z axis, describing a
distance from (or "above") the x-y plane. 
"""

torque = SymbolNew("tau", units.force * units.length, display_latex="\\tau")
"""
**Torque** is the turning effect of a force applied to a rotational system at a distance from the axis of
rotation.
"""

torsion_stiffness = SymbolNew(
    "kappa",
    units.force * units.length / angle_type,
    display_latex="\\kappa")
"""
**Torsion stiffness** or **torsion elastic modulus** is equal to the change in torque required to twist
the spring through an angle of 1 radian.

**Links:**

#. `Torsion coefficient <https://en.wikipedia.org/wiki/Torsion_spring#Torsion_coefficient>`__.
"""

bulk_modulus = SymbolNew("K", units.pressure)
"""
**Bulk modulus** of a substance is a measure of the resistance of a substance to bulk compression.
"""

poisson_ratio = SymbolNew("nu", dimensionless, display_latex="\\nu")
"""
**Poisson's ratio** is a measure of the Poisson effect, the deformation (expansion or contraction) of
a material in directions perpendicular to the specific direction of loading.
"""

engineering_normal_strain = SymbolNew("e", dimensionless)
"""
**Engineering strain**, also known as **Cauchy strain**, is expressed as the ratio of total deformation
to the initial dimension of the material body on which forces are applied.
"""

deformation = SymbolNew("Delta(l)", units.length, display_latex="\\Delta l")
"""
**Deformation** is a change in an object's shape or form due to the application of a force or forces. 
"""

strain = SymbolNew("e", dimensionless)
"""
**Strain** is defined as relative deformation, compared to a reference position configuration.
"""
