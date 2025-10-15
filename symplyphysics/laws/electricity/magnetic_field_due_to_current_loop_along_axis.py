"""
Magnetic field due to current loop along axis
=============================================

Using the Biot—Savart law, we can calculate the magnetic field due to a current loop along its
axis of symmetry. The magnetic field is directed along that axis, is proportional to the current
in the loop and depends on the distance to the center of the loop and the radius of the loop.

**Conditions:**

#. The medium is vacuum.

**Links:**

#. `Physics LibreTexts — Magnetic Field of a Current Loop <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/12%3A_Sources_of_Magnetic_Fields/12.05%3A_Magnetic_Field_of_a_Current_Loop#Eq.+12.15>`__.

..
    TODO: link to Biot-Savart law
"""

from sympy import (Eq, Rational)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities)

magnetic_flux_density = symbols.magnetic_flux_density
"""
Radial component (i.e. along the axis of the loop) of the :symbols:`magnetic_flux_density` vector.
"""

current = symbols.current
"""
Electric :symbols:`current` in the loop.
"""

loop_radius = symbols.radius
"""
:symbols:`radius` of the loop.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` from the point at which the magnetic field is measured to the center
of the loop.
"""

law = Eq(
    magnetic_flux_density,
    (quantities.vacuum_permeability * current * loop_radius**2) / (2 *
    (distance**2 + loop_radius**2)**Rational(3, 2)),
)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from Biot—Savart law


@validate_input(
    current_=current,
    loop_radius_=loop_radius,
    distance_=distance,
)
@validate_output(magnetic_flux_density)
def calculate_magnetic_flux_density(
    current_: Quantity,
    loop_radius_: Quantity,
    distance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        current: current_,
        loop_radius: loop_radius_,
        distance: distance_,
    })

    return Quantity(result)


# UNIQUE_LAW_ID: 477
