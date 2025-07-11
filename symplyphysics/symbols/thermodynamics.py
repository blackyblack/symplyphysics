"""
Thermodynamics (Symbols)
========================

Symbols related to thermodynamics.
"""

from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import Symbol

temperature = Symbol("T", units.temperature)
"""
Temperature is a scalar quantity that quantitatively expresses the attribute of hotness and coldness. It reflects
the average kinetic energy of the vibrating and colliding atoms making up a substance.
"""

adiabatic_index = Symbol("gamma", dimensionless, display_latex="\\gamma")
"""
**Adiabatic index**, or **heat capacity ratio**, is the ratio of heat capacity at constant pressure to that
at constant volume.
"""

heat_capacity = Symbol("C", units.energy / units.temperature)
"""
**Heat capacity** or **thermal capacity** is a physical property of matter, defined as the amount of heat to be
supplied to an object to produce a unit change in its temperature.
"""

molar_heat_capacity = Symbol("c_m", units.energy / (units.temperature * units.amount_of_substance))
"""
**Molar heat capacity** is defined as the heat capacity per unit amount of substance.
"""

thermal_expansion_coefficient = Symbol("alpha", 1 / units.temperature, display_latex="\\alpha")
"""
**Thermal expansion coefficient** describes how the size of an object changes with a change in temperature at
constant pressure.
"""

thermodynamic_compressibility = Symbol("beta", 1 / units.pressure, display_latex="\\beta")
"""
**Compressibility** is a measure of the instantaneous relative volume change of a fluid or solid as a response
to a pressure or mean stress change.
"""

thermal_resistance = Symbol("R", units.temperature / units.power)
"""
**Thermal resistance** measures the opposition to the heat current in a material or system.
"""

thermal_conductivity = Symbol("k", units.power / (units.length * units.temperature))
"""
**Thermal conductivity** of a material is a measure of its ability to conduct heat. It is defined as the
proportionality coefficient between the heat flux and the temperature gradient.
"""

thermal_insulance = Symbol("R_val",
    units.area * units.temperature / units.power,
    display_latex="R_\\text{val}")
"""
The **R-value**, or **thermal insulance**, is a measure of how well a two-dimensional barrier, such as a layer of
insulation, a window or a complete wall or ceiling, resists the conductive flow of heat, in the context of construction.
"""

compressibility_factor = Symbol("Z", dimensionless)
"""
The **compressibility factor**, also known as the **compression factor** or the **gas deviation factor**,
describes the deviation of a real gas from ideal gas behavior.
"""

partition_function = Symbol("Z", dimensionless)
"""
In statistical mechanics, the **partition function** describes the statistical properties of a system in
thermodynamic equilibrium. It plays the role of a normalization constant in microstate distributions of
the system by encoding the information about how the probabilities are partitioned among the different
microstates based on the specific microstate variables.
"""

boltzmann_factor = Symbol("f", dimensionless)
"""
In statistical mechanics, the **Boltzmann factor** is a quantity that describes the approximate fraction
of particles in the canonical ensemble.
"""

entropy = Symbol("S", units.energy / units.temperature)
"""
**Entropy** is a physical quantity most commonly associated with a state of randomness or disorder. In the
approach of the classical thermodynamics, entropy is defined in terms of macroscopically measurable physical
properties, such as volume, bulk mass, pressure, etc. The statistical definition defines it in terms of the
statistics of the motions of the microscopic constituents of a system.
"""

chemical_potential = Symbol("mu", units.energy, display_latex="\\mu")
"""
The **chemical potential** of a species is the energy that can be absorbed or released due to a change of the
particle number of the given species, e.g. in a chemical reaction or phase transition.
"""

gibbs_energy = Symbol("G", units.energy)
"""
The **Gibbs energy** is a thermodynamic potential that can be used to calculate the maximum amount of work,
other than pressure-volume work, that may be performed by a thermodynamically closed system at constant
temperature and pressure.
"""

enthalpy = Symbol("H", units.energy)
"""
**Enthalpy** is a state function defined as the sum of a thermodynamic system's internal energy and the
product of its pressure and volume, used in measurements at a constant external pressure.
"""

helmholtz_free_energy = Symbol("F", units.energy)
"""
In thermodynamics, the **Helmholtz free energy** (or **Helmholtz energy**) is a thermodynamic potential
that measures the useful work obtainable from a closed thermodynamic system at a constant temperature.
"""

internal_energy = Symbol("U", units.energy)
"""
**Internal energy** is a thermodynamical state function which denotes the entire energy of a closed system
of molecules or the sum of a substance's molecular kinetic and potential energy. It excludes the potential
and kinetic energies of the system as a whole and is only concerned with the energy of the molecules comprising
the system.
"""

thermal_wavelength = Symbol("lambda", units.length, display_latex="\\lambda")
"""
The **thermal de Broglie wavelength** is a quantity that is roughly the average de Broglie wavelength of
particles in an ideal gas at the specified temperature.

**Links:**

#. `Thermal de Broglie wavelength <https://en.wikipedia.org/wiki/Thermal_de_Broglie_wavelength>`__.
"""

heat = Symbol("Q", units.energy)
"""
In thermodynamics, **heat** is energy in transfer between a thermodynamic system and its surroundings by modes
other than thermodynamic work and transfer of matter.
"""

thermal_efficiency = Symbol("eta", dimensionless, display_latex="\\eta")
"""
The **thermal efficiency** is a dimensionless performance measure of a device that uses thermal energy. A generic
definition of thermal energy is the ratio of the energy benefit to the energy costs attributed to the device.
"""

statistical_weight = Symbol("Omega", dimensionless, display_latex="\\Omega")
"""
**Statistical weight**, or **multiplicity**, is a physical quantity denoting the number of microstates
corresponding to a particular macrostate of a thermodynamic system.
"""

relative_humidity = Symbol("phi", dimensionless, display_latex="\\varphi")
"""
**Relative humidity** is a quantity that indicates a present state of absolute humidity relative to a
maximum humidity given the same temperature.
"""

mean_free_path = Symbol("lambda", units.length, display_latex="\\lambda")
"""
**Mean free path** is the average distance over which a moving particle travels before
substantially changing its direction or energy.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mean_free_path>`__.
"""

sutherland_constant = Symbol("S", units.temperature)
"""
The **Sutherland constant** is a parameter used in the Sutherland model for gaseous viscosity.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Temperature_dependence_of_viscosity#Sutherland_model>`__.
"""

prandtl_number = Symbol("Pr", dimensionless, display_latex="\\text{Pr}")
"""
The **Prandtl number** is a dimensionless number defined as the ratio of momentum
diffusivity to thermal diffusivity.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Prandtl_number>`__.
"""

grashof_number = Symbol("Gr", dimensionless, display_latex="\\text{Gr}")
"""
The **Grashof number** approximates the ratio of the buoyancy to viscous forces acting
on a fluid.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Grashof_number>`__.
"""

kinematic_viscosity = Symbol("nu", units.area / units.time, display_latex="\\nu")
"""
In fluid dynamics, **kinematic viscosity**, or **momentum diffusivity**, is defined as
the ratio of the dynamic viscosity over the density of the fluid.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Viscosity#Kinematic_viscosity>`__.
"""

thermal_diffusivity = Symbol("alpha", units.area / units.time, display_latex="\\alpha")
"""
**Thermal diffusivity** is the thermal conductivity divided by density and specific heat
capacity at constant pressure.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Thermal_diffusivity>`__.
"""

attractive_forces_parameter = Symbol("a",
    units.pressure * (units.volume / units.amount_of_substance)**2)
"""
Parameter specific to each individual substance, usually attributed to the magnitude of
attractive forces between particles of the system.
"""

excluded_volume_parameter = Symbol("b", units.volume / units.amount_of_substance)
"""
Parameter specific to each individual substance, usually attributed to the amount of
excluded molar volume due to a finite size of particles.
"""

heat_transfer_coefficient = Symbol("h", units.power / (units.area * units.temperature))
"""
The **heat transfer coefficient** is the proportionality constant between the heat flux
:math:`d \\dot{Q} / d A` and the thermodynamic driving force for the flow of heat
(i.e. the temperature difference :math:`\\Delta T`). In other words, it is the rate of
heat transfer through a surface of unit area that results in a unit temperature
difference.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Heat_transfer_coefficient#>`__.
"""
