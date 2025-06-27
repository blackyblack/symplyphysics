"""
Potential energy of electric dipole in uniform electric field
=============================================================

An electric dipole placed in an external uniform electric field possesses potential energy. The
less the angle between the dipole and the electric field is, the less is its potential energy.

**Conditions:**

#. The electric field must be spatially uniform across the small region occupied by the dipole.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electric_dipole_moment#Energy_and_torque>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

potential_energy = symbols.potential_energy
"""
:symbols;`potential_energy` of the electric dipole.
"""

electric_dipole_moment = clone_as_vector_symbol(symbols.electric_dipole_moment)
"""
Vector of the :symbols:`electric_dipole_moment`.
"""

electric_field = clone_as_vector_symbol(symbols.electric_field_strength)
"""
Vector of the electric field. See :symbols:`electric_field_strength`.
"""

law = Eq(potential_energy, -1 * VectorDot(electric_dipole_moment, electric_field))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    electric_dipole_moment_=electric_dipole_moment,
    electric_field_=electric_field,
)
@validate_output(potential_energy)
def calculate_potential_energy(
    electric_dipole_moment_: QuantityCoordinateVector,
    electric_field_: QuantityCoordinateVector,
) -> Quantity:
    result = law.rhs.subs({
        electric_dipole_moment: electric_dipole_moment_,
        electric_field: electric_field_,
    })

    return Quantity(result)
