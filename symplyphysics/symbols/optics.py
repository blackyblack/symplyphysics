"""
Optics (Symbols)
================

Symbols related to optics.
"""

from sympy.physics import units
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import Symbol

relative_refractive_index = Symbol("n", dimensionless)
"""
**Relative refractive index** of of an optical medium is a dimensionless number that gives the
indication of the light bending ability of that medium. It is defined relative to a certain medium.
"""

radiant_exitance = Symbol("M_e", units.power / units.area, display_latex="M_\\text{e}")
"""
**Radiant exitance** or **radiant emittance** is the radiant flux emitted by a surface per unit area.
"""

radiant_flux = Symbol("Phi_e", units.power, display_latex="\\Phi_\\text{e}")
"""
**Radiant flux** or **radiant power** is the radiant energy emitted, reflected, transmitted, or
received per unit time.
"""

focal_length = Symbol("f", units.length)
"""
The **focal length** of an optical system is a measure of how strongly the system converges or diverges
light; it is the inverse of the system's optical power. It is the distance between the focal plane and
the lens's nodal point.
"""

optical_distance = Symbol("Lambda", units.length, display_latex="\\Lambda")
"""
**Optical distance**, also called **optical path length**, is the length that light needs to travel
through a vacuum to create the same phase difference as it would have when traveling through a given medium.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Optical_path_length>`__.
"""

irradiance = Symbol("E_e", units.power / units.area, display_latex="E_\\text{e}")
"""
In radiometry, **irradiance** is the radiant flux received by a surface per unit area.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Irradiance>`__.
"""

transparency_coefficient = Symbol("k", dimensionless)
"""
The transparency coefficient is used to describe how much light an imperfect polarizer absorbs,
zero being a total absorption and one being a perfect polarizer that lets all light pass through.

..
    TODO check if this is the correct English name
"""

reflectance = Symbol("R", dimensionless)
"""
The **reflectance** of the surface of a material is its effectiveness in reflecting radiant energy.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Reflectance>`__.
"""

magnification = Symbol("M", dimensionless)
"""
Optical **magnification** is the ratio between the apparent size in an image and its true size.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Magnification#Size_ratio_(optical_magnification)>`__.
"""

optical_power = Symbol("D", 1 / units.length)
"""
**Optical power**, also callled **dioptric power** or **focusing power**, is the degree to which
a lens, mirror, or other optical system converges or diverges light.
"""

angular_resolution = Symbol("theta", angle_type, display_latex="\\theta")
"""
**Angular resolution** describes the ability of any image-forming device such as an optical or radio telescope, a microscope, a camera, or an eye, to distinguish small details of an object.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Angular_resolution#>`__.
"""
