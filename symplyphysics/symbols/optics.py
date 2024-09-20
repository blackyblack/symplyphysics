"""
Optics (Symbols)
================

Symbols related to optics.
"""

from sympy.physics import units
from ..core.symbols.symbols import SymbolNew
from ..core.dimensions import dimensionless

relative_refractive_index = SymbolNew("n", dimensionless)
"""
**Relative refractive index** of of an optical medium is a dimensionless number that gives the
indication of the light bending ability of that medium. It is defined relative to a certain medium.
"""
