"""
Basic (Symbols)
===============

Symbols of fundamental physical quantities.
"""

from sympy.physics import units
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from symplyphysics.core.dimensions import dimensionless, any_dimension
from symplyphysics.core.symbols.symbols import SymbolNew

any_quantity = SymbolNew("X", any_dimension)
"""
A quantity that can have any dimension.
"""

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

molar_mass = SymbolNew("M", units.mass / units.amount_of_substance)
"""
**Molar mass** is defined as the mass per unit amount of substance.
"""

work = SymbolNew("W", units.energy)
"""
**Work** is the energy transferred to or from an object via the application of force.
"""

energy_density = SymbolNew("w", units.energy / units.volume)
"""
Energy per unit volume.
"""

spectral_energy_density = SymbolNew("w_f", units.energy / (units.volume * units.frequency))
"""
Energy per unit volume per unit linear frequency.
"""

energy = SymbolNew("E", units.energy)
"""
**Energy** is the quantitative property that is transferred to a body or to a physical system, recognizable in
the performance of work and in the form of heat and light.
"""

specific_energy = SymbolNew("epsilon", units.energy / units.mass, display_latex="\\varepsilon")
"""
**Specific energy** is defined as :symbols:`energy` per unit :symbols:`mass`.
"""

power = SymbolNew("P", units.power)
"""
**Power** is the amount of energy transferred or converted per unit time.
"""

radius_of_curvature = SymbolNew("r", units.length)
"""
**Radius of curvature** is the inverse of curvature and is equal to the distance to the center of curvature.

..
    TODO are this and distance_to_axis interchangeable?
"""

density = SymbolNew("rho", units.mass / units.volume, display_latex="\\rho")
"""
**Density** is mass per unit volume.
"""

linear_density = SymbolNew("mu", units.mass / units.length, display_latex="\\mu")
"""
**Linear density** is mass per unit length.
"""

intensity = SymbolNew("I", units.power / units.area)
"""
**Intensity** or **flux** of radiant energy is the power transferred per unit area,  where the area is measured
on the plane perpendicular to the direction of propagation of the energy.
"""

whole_number = SymbolNew("N", dimensionless, integer=True)
"""
A dimensionless **whole** number of any sign.
"""

positive_number = SymbolNew("N", dimensionless, integer=True, positive=True)
"""
A dimensionless whole **number** used for counting objects or instances.
"""

nonnegative_number = SymbolNew("N", dimensionless, integer=True, nonnegative=True)
"""
A dimensionless non-negative whole **number**, i.e. :math:`0, 1, 2, \\dots`.
"""

number_density = SymbolNew("n", 1 / units.volume)
"""
**Number density** is an intensive quantity used to describe the degree of concentration of countable objects
(particles, molecules, phonons, cells, galaxies, etc.) in physical space.
"""

particle_count = SymbolNew("N", dimensionless)
"""
Number of particles in the system.
"""

angle = SymbolNew("phi", angle_type, display_latex="\\varphi")
"""
An **angle** is the difference in direction between two lines or surfaces.
"""

probability = SymbolNew("P", dimensionless)
"""
**Probability** is a measure of an event's likelihood.
"""

fractional_change = SymbolNew("e", dimensionless)
"""
**Fractional change** is linear change divided by initial value of the quantity.
"""

exponential_decay_constant = SymbolNew("lambda", 1 / units.time, display_latex="\\lambda")
"""
**Exponential decay constant**, also called **rate constant** or **disintegration constant**, is 
the rate at which some quantity is decreasing in such a way that its rate of change is proportional
to its current value.

**Links:**

#. `Exponential decay <https://en.wikipedia.org/wiki/Exponential_decay>`__.
"""

characteristic_length = SymbolNew("l_c", units.length, display_latex="l_\\text{c}")
"""
**Characteristic length** is a dimension that defines the scale of the physical system.
It is usually defined as the volume of the system divided by its surface.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Characteristic_length>`__.
"""

molar_volume = SymbolNew("v_m", units.volume / units.amount_of_substance, display_latex="v_\\text{m}")
"""
**Molar volume** is defined as :symbols:`volume` per unit :symbols:`amount_of_substance`.
"""
