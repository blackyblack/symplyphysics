"""
Chemistry (Symbols)
===================

Symbols of physical quantities related to chemistry.
"""

from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import SymbolNew

mass_fraction = SymbolNew("w", dimensionless)
"""
**Mass fraction** of a substance within a mixture is the ratio of the mass of the substance to the total mass
of the mixture.
"""

amount_of_substance = SymbolNew("n", units.amount_of_substance)
"""
**Amount of substance** in a given sample of matter is defined as a ratio between the number of elementary
entities and the Avogadro constant.
"""

density_of_states = SymbolNew("D", 1 / units.volume)
r"""
The **density of states** of a system describes describes the number of allowed modes or states
per unit energy range.
"""

band_gap = SymbolNew("E_g", units.energy, display_latex="E_\\text{g}")
"""
A **band gap**, or **energy gap**, is an energy range in a solid where no electronic states exist. 
"""

work_function = SymbolNew("W", units.energy)
"""
**Work function** is the minimum thermodynamic work (i.e., energy) needed to remove an electron from
a solid to a point in the vacuum immediately outside the solid surface.
"""

drift_velocity = SymbolNew("u", units.velocity)
"""
**Drift velocity** is the average velocity attained by charged particles, such as electrons, in a material
due to an electric field.
"""

molar_concentration = SymbolNew("c", units.amount_of_substance / units.volume)
"""
**Molar concentration**, or **molarity**, is a quantity most commonly defined as amount of substance of
solute per unit volume of solution, or per unit volume available to the species.
"""

ionization_coefficient = SymbolNew("alpha", 1 / units.length, display_latex="\\alpha")
"""
**Ionization coefficient** can be defined as the mean number of ionization processes
over the distance covered in the direction of the electric field.

**Links:**

#. `ETH Research Collection <https://www.research-collection.ethz.ch/bitstream/20.500.11850/186582/1/PostPrint.pdf>`__.
"""

cross_section = SymbolNew("sigma", units.area, display_latex="\\sigma")
"""
**Cross section** is a measure of the probability that a specific process will take place
in a collision of two particles.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Cross_section_(physics)>`__.
"""
