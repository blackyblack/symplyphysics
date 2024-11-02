"""
Thermodynamics (Symbols)
========================

Symbols related to thermodynamics.
"""

from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import SymbolNew

temperature = SymbolNew("T", units.temperature)
"""
Temperature is a scalar quantity that quantitatively expresses the attribute of hotness and coldness. It reflects
the average kinetic energy of the vibrating and colliding atoms making up a substance.
"""

adiabatic_index = SymbolNew("gamma", dimensionless, display_latex="\\gamma")
"""
**Adiabatic index**, or **heat capacity ratio**, is the ratio of heat capacity at constant pressure to that
at constant volume.
"""

heat_capacity = SymbolNew("C", units.energy / units.temperature)
"""
**Heat capacity** or **thermal capacity** is a physical property of matter, defined as the amount of heat to be
supplied to an object to produce a unit change in its temperature.
"""

molar_heat_capacity = SymbolNew("c_m", units.energy / (units.temperature * units.amount_of_substance))
"""
**Molar heat capacity** is defined as the heat capacity per unit amount of substance.
"""

thermal_expansion_coefficient = SymbolNew("alpha", 1 / units.temperature, display_latex="\\alpha")
"""
**Thermal expansion coefficient** describes how the size of an object changes with a change in temperature at
constant pressure.
"""

thermodynamic_compressibility = SymbolNew("beta", 1 / units.pressure, display_latex="\\beta")
"""
**Compressibility** is a measure of the instantaneous relative volume change of a fluid or solid as a response
to a pressure or mean stress change.
"""

thermal_resistance = SymbolNew("R", units.temperature / units.power)
"""
**Thermal resistance** measures the opposition to the heat current in a material or system.
"""

thermal_conductivity = SymbolNew("k", units.power / (units.length * units.temperature))
"""
**Thermal conductivity** of a material is a measure of its ability to conduct heat. It is defined as the
proportionality coefficient between the heat flux and the temperature gradient.
"""

thermal_insulance = SymbolNew("R_val",
    units.area * units.temperature / units.power,
    display_latex="R_\\text{val}")
"""
The **R-value**, or **thermal insulance**, is a measure of how well a two-dimensional barrier, such as a layer of
insulation, a window or a complete wall or ceiling, resists the conductive flow of heat, in the context of construction.
"""

compressibility_factor = SymbolNew("Z", dimensionless)
"""
The **compressibility factor**, also known as the **compression factor** or the **gas deviation factor**,
describes the deviation of a real gas from ideal gas behavior.
"""

partition_function = SymbolNew("Z", dimensionless)
"""
In statistical mechanics, the **partition function** describes the statistical properties of a system in
thermodynamic equilibrium. It plays the role of a normalization constant in microstate distributions of
the system by encoding the information about how the probabilities are partitioned among the different
microstates based on the specific microstate variables.
"""

boltzmann_factor = SymbolNew("f", dimensionless)
"""
In statistical mechanics, the **Boltzmann factor** is a quantity that describes the approximate fraction
of particles in the canonical ensemble.
"""

entropy = SymbolNew("S", units.energy / units.temperature)
"""
Entropy is a physical quantity most commonly associated with a state of randomness or disorder. In the
approach of the classical thermodynamics, entropy is defined in terms of macroscopically measurable physical
properties, such as volume, bulk mass, pressure, etc. The statistical definition defines it in terms of the
statistics of the motions of the microscopic constituents of a system.
"""
