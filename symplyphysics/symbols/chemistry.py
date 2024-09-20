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
