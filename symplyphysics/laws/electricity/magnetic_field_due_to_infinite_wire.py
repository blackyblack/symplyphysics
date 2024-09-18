r"""
Magnetic field due to infinite wire
===================================

The magnitude of the magnetic flux density due to a thin, straight, infinite wire depends on the
current through it and the radial distance to the wire.

**Notation:**

#. :math:`\mu_0` (:code:`mu_0`) is vacuum permeability.

**Conditions:**

#. The wire is uniform, straight, and thin.
"""

from sympy import (Eq, pi)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

magnetic_flux_density = symbols.magnetic_flux_density
"""
Magnitude of :symbols:`magnetic_flux_density`.
"""

current = symbols.current
"""
:symbols:`current` flowing through the wire.
"""

radial_distance = symbols.distance_to_axis
"""
Radial distance to wire. See :symbols:`distance_to_axis`.
"""

law = Eq(magnetic_flux_density, quantities.vacuum_permeability * current / (2 * pi * radial_distance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    current_=current,
    radial_distance_=radial_distance,
)
@validate_output(magnetic_flux_density)
def calculate_magnetic_field(
    current_: Quantity,
    radial_distance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        current: current_,
        radial_distance: radial_distance_,
    })
    return Quantity(result)
