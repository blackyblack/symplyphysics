"""
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
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

magnetic_flux_density = Symbol("magnetic_flux_density", units.magnetic_flux_density)
"""
Magnitude of magnetic flux density.

Symbol:
    :code:`B`
"""

current = Symbol("current", units.current)
"""
Current flowing through the wire.

Symbol:
    :code:`I`
"""

radial_distance = Symbol("radial_distance", units.length)
"""
Radial distance to wire.

Symbol:
    :code:`r`
"""

law = Eq(magnetic_flux_density, units.vacuum_permeability * current / (2 * pi * radial_distance))
r"""
:code:`B = mu_0 * I / (2 * pi * r)`

Latex:
    .. math::
        B = \frac{\mu_0 I}{2 \pi r}
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
