"""
Symbols (Basic)
===============

Contains fundamental properties of matter.
"""

from sympy.physics import units
from ...core.symbols.symbols import Symbol

mass = Symbol("mass", units.mass)
"""
Mass can be experimentally defined as a measure of the body's inertia, meaning the resistance to acceleration (change of velocity) when a net force is applied.
The object's mass also determines the strength of its gravitational attraction to other bodies.

Symbols:
    m, M
"""

__all__ = [
    "mass",
]
