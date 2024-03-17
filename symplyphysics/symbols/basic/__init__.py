"""! @brief Contains fundamental properties of matter."""

from sympy.physics import units
from ...core.symbols.symbols import Symbol

# Mass can be experimentally defined as a measure of the body's inertia, meaning the resistance to acceleration (change of velocity) when a net force is applied.
# The object's mass also determines the strength of its gravitational attraction to other bodies.
mass = Symbol("mass", units.mass)

__all__ = [
    "mass",
]
