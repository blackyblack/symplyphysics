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

radiant_exitance = SymbolNew("M_e", units.power / units.area, display_latex="M_\\text{e}")
"""
**Radiant exitance** or **radiant emittance** is the radiant flux emitted by a surface per unit area.
"""

radiant_flux = SymbolNew("Phi_e", units.power, display_latex="\\Phi_\\text{e}")
"""
**Radiant flux** or **radiant power** is the radiant energy emitted, reflected, transmitted, or
received per unit time.
"""
