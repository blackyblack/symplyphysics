"""
Astronomy (Symbols)
===================

Symbols related to astronomy.
"""

from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import SymbolNew

absolute_magnitude = SymbolNew("M", dimensionless)
"""
**Absolute magnitude** is a measure of the luminosity of a celestial object on an inverse
logarithmic astronomical magnitude scale.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Absolute_magnitude>`__.
"""

apparent_magnitude = SymbolNew("m", dimensionless)
"""
**Apparent magnitude** is a measure of the brightness of a star, astronomical object or other
celestial objects like artificial satellites. Its value depends on its intrinsic luminosity,
its distance, and any extinction of the object's light caused by interstellar dust along the
line of sight to the observer.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Apparent_magnitude>`__.
"""

luminocity = SymbolNew("L", units.power)
"""
In astronomy, **luminosity** is the total amount of electromagnetic energy emitted per unit of
time by a star, galaxy, or other astronomical objects.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Luminosity#>`__.
"""
