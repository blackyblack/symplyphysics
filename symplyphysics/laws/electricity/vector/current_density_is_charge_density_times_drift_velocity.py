"""
Current density is charge density times drift velocity
======================================================

The vector of the electric current density is a vector whose magnitude is the electric current per
cross-sectional area at a given point in space and its direction is that of the motion of positive
electric charges at that point.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Current_density>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol

current_density = clone_as_vector_symbol(symbols.current_density)
"""
Vector of the :symbols:`current_density`.
"""

charge_density = symbols.volumetric_charge_density
"""
:symbols:`volumetric_charge_density`.
"""

drift_velocity = clone_as_vector_symbol(symbols.drift_velocity)
"""
Vector of the :symbols:`drift_velocity`.
"""

law = Eq(current_density, charge_density * drift_velocity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    charge_density_=charge_density,
    drift_velocity_=drift_velocity,
)
@validate_output(current_density)
def calculate_current_density(
    charge_density_: Quantity,
    drift_velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        charge_density: charge_density_,
        drift_velocity: drift_velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 530
