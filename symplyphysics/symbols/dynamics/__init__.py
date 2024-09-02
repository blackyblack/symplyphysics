"""
Symbols (Dynamics)
==================

Contains properties for 'dynamics' branch of physics.

Dynamics is subdivision of mechanics that is concerned with the motion of material objects in relation to the physical factors that affect them.
"""

from sympy.physics import units
from ...core.symbols.symbols import SymbolNew

force = SymbolNew("F", units.force)
"""
Force is an influence that can cause an object to change its velocity, i.e., to accelerate, meaning a change in speed or direction, unless counterbalanced by other forces.
"""

__all__ = [
    "force",
]
