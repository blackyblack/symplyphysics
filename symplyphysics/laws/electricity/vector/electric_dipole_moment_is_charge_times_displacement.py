"""
Electric dipole moment is charge times displacement
===================================================

The vector of electric dipole moment is a vector whose magnitude describes the
:doc:`electric dipole moment <laws.electricity.electric_dipole_moment_is_charge_times_distance>`
of the system. It is collinear to the vector connecting the two point charges.

**Notes:**

#. This law can be seen as the definition of electric dipole moment.

**Conditions:**

#. The charges must be equal by magnitude and have opposite signs.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electric_dipole_moment#Elementary_definition>`__.
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, Quantity, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

electric_dipole_moment = clone_as_vector_symbol(symbols.electric_dipole_moment)
"""
Vector of :symbols:`electric_dipole_moment`.
"""

charge = symbols.charge
"""
Magnitude of the electric :symbols:`charge` of the point charges.
"""

displacement = clone_as_vector_symbol(
    symbols.distance,
    display_symbol="d",
    display_latex="{\\vec d}",
)
"""
Position vector pointing from the *negative* charge to the *positive* charge. See
:symbols:`distance`.
"""

law = Eq(electric_dipole_moment, charge * displacement)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_=charge, displacement_vector_=displacement)
@validate_output(electric_dipole_moment)
def calculate_dipole_moment(
    charge_: Quantity,
    displacement_vector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        charge: charge_,
        displacement: displacement_vector_,
    })

    return QuantityCoordinateVector.from_expr(result)
