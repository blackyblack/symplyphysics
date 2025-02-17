"""
Classical mechanics (Symbols)
=============================

Symbols related to classical mechanics.
"""

from sympy.physics import units
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import Symbol

force = Symbol("F", units.force)
"""
**Force** is a vector quantity denoting an influence that can cause an object to change its velocity
unless counterbalanced by other forces.
"""

tension = Symbol("T", units.force)
"""
**Tension** is the pulling or stretching force transmitted axially along an object so as to stretch it
or pull it apart.
"""

speed = Symbol("v", units.velocity)
"""
**Speed** is the rate of change of the object's position with respect to time. It is the scalar counterpart
of the velocity vector, which is the rate of change of the object's displacement with respect to time.
"""

acceleration = Symbol("a", units.acceleration)
"""
**Acceleration** is the rate of change of the object's velocity with respect to time. It is a vector quantity.
"""

position = Symbol("x", units.length)
"""
**Position** is defined as the spatial location of an object with respect to a coordinate system.
"""

euclidean_distance = Symbol("d", units.length)
r"""
**Euclidean distance** between two points is the length of the line segment between them. Mathematically,
it can be represented by the following formula:

.. math::
    d = \sqrt{\left\Vert \mathbf{r}_2 - \mathbf{r}_1 \right\Vert},

where :math:`\mathbf{r}_1` is the position vector of the first point, :math:`\mathbf{r}_2` is the position
vector of the second point, and :math:`\left\Vert \cdot \right\Vert` is the Euclidean norm.

**Notes:**

#. This symbol can be used for representing **displacement**. In that case, :math:`\mathbf{r}_1` corresponds
   to the initial position of the body and :math:`\mathbf{r}_2` to the final position of the body.
"""

distance = Symbol("s", units.length)
r"""
**distance**, or more precisely **total distance traveled**, is the length of the path traced by a moving body.
Mathematically, it can be represented by the following formula:

.. math::
    s(t) = \int \limits_0^t |\mathbf{v}(\tau)| d\tau,

where :math:`\mathbf{v}` is the velocity vector of the body as a function of time :math:`t`. It is assumed that
the total distance traveled is zero at zero time.

**Notes:**

#. It is sometimes used in the sense of *the length of space between two points*, in which case
   :symbols:`euclidean_distance` must be used.
"""

distance_to_origin = Symbol("r", units.length)
"""
Distance to the origin of the coordinate system.
"""

distance_to_axis = Symbol("r", units.length)
"""
Distance to reference axis.
"""

orthogonal_distance = Symbol("z", units.length)
"""
**Orthogonal distance** between two objects is the distance from one to the other,
measured along a line that is perpendicular to one or both.
"""

length = Symbol("l", units.length)
"""
**Length** is a measure of a size of an object.
"""

radius = Symbol("r", units.length)
"""
**Radius** of a sphere is the distance from the center of the sphere to any point on
the sphere.
"""

diameter = Symbol("D", units.length)
"""
**Diameter** of a circle or sphere is the length of the biggest chord connecting two
points on the circumference.
"""

semimajor_axis = Symbol("a", units.length)
"""
**Semi-major axis**, or **major semiaxis**, is the longest semidiameter of an ellipse.
"""

semiminor_axis = Symbol("b", units.length)
"""
**Semi-minor axis**, or **minor semiaxis**, is the smallest semidiameter of an ellipse.
"""

thickness = Symbol("h", units.length)
"""
**Thickness** is a measure of a size of an object, usually the separation between two layers, or
the distance through an object distinct from length and width.
"""

area = Symbol("A", units.area)
"""
**Area** is the size of a region on a two-dimensional surface.
"""

arc_length = Symbol("s", units.length)
"""
**Arc length** is the distance between two points along a section of a curve.
"""

angular_speed = Symbol("w", angle_type / units.time, display_latex="\\omega")
"""
**Angular speed** is the rate of change of angular distance with respect to time.
"""

angular_frequency = angular_speed
"""
**Angular frequency**, also called **angular rate**, is a scalar measure of the temporal
rate of change of the phase argument of a sinusoidal waveform or sine function.
"""

angular_acceleration = Symbol("alpha", angle_type / units.time**2, display_latex="\\alpha")
"""
**Angular acceleration** is the rate of change of angular speed with respect to time.
"""

angular_distance = Symbol("theta", angle_type, display_latex="\\theta")
"""
**Angular distance** is a measure of an angular separation between two points.
"""

angular_wavenumber = Symbol("k", angle_type / units.length)
"""
**Angular wavenumber** is the spatial analog of temporal frequency equal to radians per unit length.
"""

wavelength = Symbol("lambda", units.length, display_latex="\\lambda")
"""
**Wavelength** or **spatial period** is the distance over which the wave's shape repeats.
"""

damping_ratio = Symbol("zeta", dimensionless, display_latex="\\zeta")
"""
**Damping ratio** is a dimensionless measure describing how oscillations in a system decay after a disturbance.
"""

volume = Symbol("V", units.volume)
"""
Volume is a measure of regions in three-dimensional space.
"""

impulse = Symbol("J", units.momentum)
"""
**Impulse** is the change in momentum of an object.
"""

phase_speed = Symbol("v", units.speed)
"""
**Phase speed** is the speed at which the phase of the wave travels.
"""

group_speed = Symbol("v_g", units.speed, display_latex="v_\\text{g}")
"""
**Group speed** of a wave is the speed with which the overall envelope shape of the wave's
amplitudes, or the modulation or envelope of the wave, propagates through space.
"""

pressure = Symbol("p", units.pressure)
"""
**Pressure** is the force applied perpendicular to the surface of an object per unit area over which
that force is distributed.
"""

temporal_frequency = Symbol("f", units.frequency)
"""
**Temporal frequency** is the number of occurrences of a repeating event per unit of time.
"""

sound_intensity_level = Symbol("L_I", dimensionless)
"""
**Sound intensity level** is the measure of the *intensity* of a sound relative to a reference value.
"""

rotational_inertia = Symbol("I", units.mass * units.length**2)
"""
**Rotational inertia**, also known as **moment of inertia**, is defined relative to a rotational axis
and is the ratio between the torque applied and the resulting angular acceleration about that axis.
"""

quality_factor = Symbol("Q", dimensionless)
"""
**Quality factor** or **Q factor** is a dimensionless parameter that describes how underdamped an
oscillator or resonator is.
"""

momentum = Symbol("p", units.momentum)
"""
**Momentum**, more specifically **linear** or **transitional momentum**, is the product of the mass and
velocity of an object.
"""

mechanical_energy = Symbol("E", units.energy)
"""
**Mechanical energy** is defined to be the sum of potential energy and kinetic energy.
"""

kinetic_energy = Symbol("K", units.energy)
"""
**Kinetic energy** of an object is the form of energy that it possesses due to its motion.
"""

potential_energy = Symbol("U", units.energy)
"""
**Potential energy** is the energy held by an object because of its position relative to other objects,
stresses within itself, its electric charge, or other factors. Potential energy is associated with so
called conservative forces and only depends on the initial and final positions of the body in space.
"""

mass_flow_rate = Symbol("mu", units.mass / units.time, display_latex="\\mu")
"""
**Mass flow rate** is the mass of a substance which passes per unit time.
"""

volumetric_flow_rate = Symbol("Q", units.volume / units.time)
"""
**Volumetric flow rate** is the volume of a substance which passes per unit time.
"""

stiffness = Symbol("k", units.force / units.length)
"""
**Stiffness** is the extent to which an object resists deformation in response to an applied force.
"""

compliance = Symbol("c", units.length / units.force)
"""
**Compliance** is the inverse of :symbols:`stiffness`.
"""

young_modulus = Symbol("E", units.pressure)
"""
**Young modulus** is a mechanical property of solid materials that measures the tensile or compressive
stiffness when the force is applied lengthwise.
"""

phase_shift = Symbol("phi", angle_type, display_latex="\\varphi")
"""
**Phase shift**, also known as **phase offset** or **phase difference**, is the shift of phase between
two periodic functions.
"""

phase = Symbol("phi", angle_type, display_latex="\\varphi")
"""
**Phase** of a wave or other periodic function of some real variable :math:`t`, such as time, is an
angle-like quantity representing the fraction of the cycle covered up to :math:`t`.
"""

friction_coefficient = Symbol("mu", dimensionless, display_latex="\\mu")
"""
**Coefficient of friction** is a dimensionless scalar value which equals to the ratio of the force of
friction between two bodies and the force pressing them together, either during or at the onset of
slipping.
"""

height = Symbol("h", units.length)
"""
**Height** is measure of vertical distance, either vertical extent or vertical position.

In the case of three-dimensional space, height is measured along the vertical z axis, describing a
distance from (or "above") the x-y plane. 
"""

torque = Symbol("tau", units.force * units.length, display_latex="\\tau")
"""
**Torque** is the turning effect of a force applied to a rotational system at a distance from the axis of
rotation.
"""

torsion_stiffness = Symbol("kappa",
    units.force * units.length / angle_type,
    display_latex="\\kappa")
"""
**Torsion stiffness** or **torsion elastic modulus** is equal to the change in torque required to twist
the spring through an angle of 1 radian.

**Links:**

#. `Torsion coefficient <https://en.wikipedia.org/wiki/Torsion_spring#Torsion_coefficient>`__.
"""

bulk_modulus = Symbol("K", units.pressure)
"""
**Bulk modulus** of a substance is a measure of the resistance of a substance to bulk compression.
"""

poisson_ratio = Symbol("nu", dimensionless, display_latex="\\nu")
"""
**Poisson's ratio** is a measure of the Poisson effect, the deformation (expansion or contraction) of
a material in directions perpendicular to the specific direction of loading.
"""

engineering_normal_strain = Symbol("e", dimensionless)
"""
**Engineering strain**, also known as **Cauchy strain**, is expressed as the ratio of total deformation
to the initial dimension of the material body on which forces are applied.
"""

deformation = Symbol("Delta(l)", units.length, display_latex="\\Delta l")
"""
**Deformation** is a change in an object's shape or form due to the application of a force or forces. 
"""

strain = Symbol("e", dimensionless)
"""
**Strain** is defined as relative deformation, compared to a reference position configuration.
"""

stress = Symbol("sigma", units.pressure, display_latex="\\sigma")
"""
**Stress** is a physical quantity that describes forces present during deformation.
"""

mach_number = Symbol("M", dimensionless, display_latex="\\text{M}", positive=True)
"""
The **Mach number** is a dimensionless quantity in fluid dynamics representing the ratio of flow velocity
past a boundary to the local speed of sound.
"""

diffusion_coefficient = Symbol("D", units.area / units.time)
"""
**Diffusion coefficient**, or **diffusivity**, is usually defined as the proportionality constant between
the molar flux due to molecular diffusion and the negative value of the gradient in the concentration of
the species.

**Links:**

#. `Mass diffusivity <https://en.wikipedia.org/wiki/Mass_diffusivity>`__.
"""

dynamic_viscosity = Symbol("mu", units.pressure * units.time)
"""
**Dynamic viscosity** is a physical quantity measuring the resistance to deformation at a given rate.
Specifically in fluid mechanics, it is the proportionality factor between the shear stress of the adjacent
layers of the fluid and the local sheer velocity.

**Links:**

#. `Dynamic viscosity <https://en.wikipedia.org/wiki/Viscosity#Dynamic_viscosity>`__.
"""

diffusion_flux = Symbol("J", units.amount_of_substance / (units.area * units.time))
"""
**Diffusion flux** is a physical quantity that measures the amount of substance that will flow through
a unit area during a unit time interval. For the general definition of flux, see `Flux
<https://en.wikipedia.org/wiki/Flux>`__.
"""

degrees_of_freedom = Symbol("f", dimensionless, integer=True, positive=True)
"""
A **degree of freedom** is a physical parameter in the parameterization of a physical system. The number
of degrees of freedom indicates the smallest number of parameters whose values determine all parameters
in the chosen parameterization.
"""

angular_momentum = Symbol("L", units.mass * units.length**2 / units.time)
"""
**Angular momentum**, sometimes called **rotational momentum**, is the rotational analog of linear
:symbols:`momentum`.
"""

latitude = Symbol("phi", angle_type, display_latex="\\phi")
"""
**Latitude** is a coordinate that specifies the north-south position of a point on the surface of
the Earth or another celestial body. Its value ranges from :math:`-90^\\circ` at the south pole
to :math:`90^\\circ` at the north pole, with :math:`0^\\circ` at the Equator.
"""

longitude = Symbol("lambda", angle_type, display_latex="\\lambda")
"""
**Longitude** is a coordinate that specifies the east-west position of a point on the
surface of the Earth, or another celestial body. The prime meridian defines :math:`0^\\circ`
longitude; positive longitudes are east of the prime meridian, and the negative ones are west.
"""

sector_speed = Symbol("sigma", units.area / units.time, display_latex="\\sigma")
"""
**Areal speed**, also called **sector speed** or **sectorial speed**, is a quantity that indicates
the rate of change at which :symbols:`area` is swept out by a particle as it moves along a curve.
The calculation of the area swept can be found in the link below.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Areal_velocity#/media/File:ArealVelocity_with_curved_area.svg>`__.
"""

kepler_constant = Symbol("K", units.length**3 / units.time**2, display_latex="\\mathfrak{K}")
"""
The **Kepler's constant** is the constant of proportionality in Kepler's third law of planetary
motion, namely the ratio between the square of the period of the planet to the semi-major axis
of the planet's orbit. It is constant for all objects orbiting around the same object.
"""

surface_tension = Symbol("gamma", units.force / units.length, display_latex="\\gamma")
"""
**Surface tension**, as a phenomenon, is the tendency of liquid surfaces at rest to
shrink into the minimum surface area possible. The coefficient of surface tension is
directly proportional to the force necessary to increase the surface area of the fluid
and inversely proportional to the length of the side at which the force is applied to.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Surface_tension#Physics>`__.
"""

dynamic_pressure = Symbol("q", units.pressure)
"""
**Dynamic pressure** is the :symbols:`kinetic_energy` per unit :symbols:`volume` of a
fluid. See :ref:`Quantity is volumetric density times volume`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Dynamic_pressure#Physical_meaning>`__.
"""

mechanical_efficiency = Symbol("eta", dimensionless, display_latex="\\eta")
"""
**Mechanical efficiency** is a dimensionless ratio that measures the efficiency of a
mechanism or machine in transforming the power input to the device to power output.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mechanical_efficiency#>`__.
"""

hydrostatic_pressure = Symbol("p", units.pressure)
"""
**Hydrostatic pressure** is the pressure exerted by a fluid at equilibrium at a given
point within the fluid, due to the force of gravity.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Hydrostatics#Hydrostatic_pressure>`__.
"""

froude_number = Symbol("Fr", dimensionless, display_latex="\\text{Fr}")
"""
The **Froude number** is a dimensionless number defined as the ratio of the flow inertia
to the external force field.
"""

speed_of_sound = Symbol("c", units.velocity)
"""
The **speed of sound** is the distance traveled per unit of time by a sound wave as it
propagates through an elastic medium.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Speed_of_sound#>`__.
"""

nusselt_number = Symbol("Nu", dimensionless, display_latex="\\text{Nu}")
"""
The **Nusselt number** is the ratio of total heat transfer (convection + conduction) to
conductive heat transfer across a boundary of a fluid.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Nusselt_number>`__.
"""

reynolds_number = Symbol("Re", dimensionless, display_latex="\\text{Re}")
"""
The **Reynolds number** is a dimensionless quantity that helps predict fluid flow
patterns in different situations by measuring the ratio between inertial and viscous
forces.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Reynolds_number#>`__.
"""

shear_stress = Symbol("tau", units.pressure, display_latex="\\tau")
"""
**Shear stress** is the component of stress coplanar with a material cross section.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Shear_stress>`__.
"""

flow_speed = Symbol("u", units.velocity)
"""
**Flow speed** is the magnitude of the vector of **flow velocity** (also called
**macroscopic velocity**). Flow velocity is a vector field used to mathematically
describe the motion of a continuum.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Flow_velocity>`__.
"""
