"""
Relativistic mechanics (Symbols)
================================

Symbols related to relativistic mechanics.
"""

from sympy.physics import units
from ..core.symbols.symbols import SymbolNew
from ..core.dimensions import dimensionless

lorentz_factor = SymbolNew("gamma", dimensionless, display_latex="\\gamma")
"""
**Lorentz factor** is a quantity expressing how much the measurements of time, length, and other physical
properties change for an object while it moves.
"""
