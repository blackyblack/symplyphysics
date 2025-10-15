"""
Torque due to electric dipole moment in uniform electric field
==============================================================

An object with an electric dipole moment is subject to a torque when placed in an external
electric field. The torque tends to align the dipole with the field.

**Notes:**

#. The dipole receives *no* overall net force due to the fact that the forces exerted on the point
   charges constituting the dipole cancel each other.

**Conditions:**

#. The electric field must be spatially uniform across the small region occupied by the dipole.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electric_dipole_moment#Energy_and_torque>`__.
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

torque = clone_as_vector_symbol(symbols.torque)
"""
Pseudovector of the :symbols:`torque` acting on the object.
"""

electric_dipole_moment = clone_as_vector_symbol(symbols.electric_dipole_moment)
"""
Vector of the :symbols:`electric_dipole_moment` of the object.
"""

electric_field = clone_as_vector_symbol(symbols.electric_field_strength)
"""
Vector of the uniform electric field. See :symbols:`electric_field_strength`.
"""

law = Eq(torque, VectorCross(electric_dipole_moment, electric_field))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    electric_dipole_moment_=electric_dipole_moment,
    electric_field_=electric_field,
)
@validate_output(torque)
def calculate_torque(
    electric_dipole_moment_: QuantityCoordinateVector,
    electric_field_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        electric_dipole_moment: electric_dipole_moment_,
        electric_field: electric_field_,
    })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 534
