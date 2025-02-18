"""
Astronomy (Symbols)
===================

Symbols related to astronomy.
"""

from sympy.physics import units
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import Symbol

absolute_magnitude = Symbol("M", dimensionless)
"""
**Absolute magnitude** is a measure of the luminosity of a celestial object on an inverse
logarithmic astronomical magnitude scale. An object's absolute magnitude is defined to be
equal to the apparent magnitude that the object would have if it were viewed from a
distance of exactly 10 parsecs, without extinction or dimming of its light due to
absorption by interstellar matter and cosmic dust.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Absolute_magnitude>`__.
"""

apparent_magnitude = Symbol("m", dimensionless)
"""
**Apparent magnitude** is a measure of the brightness of a star, astronomical object or other
celestial objects like artificial satellites. Its value depends on its intrinsic luminosity,
its distance, and any extinction of the object's light caused by interstellar dust along the
line of sight to the observer.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Apparent_magnitude>`__.
"""

luminosity = Symbol("L", units.power)
"""
In astronomy, **luminosity** is the total amount of electromagnetic energy emitted per unit of
time by a star, galaxy, or other astronomical objects.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Luminosity#>`__.
"""

illuminance = Symbol("E_v", units.luminous_intensity / units.area, display_latex="E_\\text{v}")
"""
In photometry, **illuminance** is the total luminous flux incident on a surface, per unit area.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Illuminance#>`__.
"""

irradiance = Symbol("E_e", units.power / units.area, display_latex="E_\\text{e}")
"""
In radiometry, **irradiance**, or **flux density**, is the radiant flux received by a
surface per unit area.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Irradiance>`__.
"""

eccentricity = Symbol("e", dimensionless)
"""
The **eccentricity** of a conic section is a non-negative real number that uniquely
characterizes its shape. 
"""

zenith_angle = Symbol("theta", angle_type)
"""
The **zenith angle**, or **zenith angular distance**, is the angle between a direction of
interest (e.g. a star) and the local zenith.
"""

declination = Symbol("delta", angle_type, display_latex="\\delta")
"""
**Declination** is one of the two angles that locate a point on the celestial sphere in
the equatorial coordinate system, measured north (positive) or south (negative) of the
celestial equator.
"""

altitude = Symbol("h", angle_type)
"""
**Altitude**, sometimes referred to as **elevation** or **apparent height**, is the angle
between the object and the observer's local horizon.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Horizontal_coordinate_system#Definition>`__.
"""
