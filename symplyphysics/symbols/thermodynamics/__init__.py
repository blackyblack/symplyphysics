"""
Symbols (Thermodynamics)
========================

Contains properties for 'thermodynamics' branch of physics

Thermodynamics is a branch of physics that deals with heat, work, and temperature, and their relation to energy, entropy, and the physical properties of matter and radiation.
"""

from sympy.physics import units
from ...core.symbols.symbols import SymbolNew

temperature = SymbolNew("T", units.temperature)
"""
Temperature reflects the average kinetic energy of the vibrating and colliding atoms making up a substance.
"""

__all__ = [
    "temperature",
]
