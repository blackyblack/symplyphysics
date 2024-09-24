"""
Chemistry (Symbols)
===================

Symbols of physical quantities related to chemistry.
"""

from sympy.physics import units
from ..core.dimensions import dimensionless
from ..core.symbols.symbols import SymbolNew

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
