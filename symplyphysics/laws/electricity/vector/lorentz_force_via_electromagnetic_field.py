"""
Lorentz force via electromagnetic field
=======================================

The **Lorentz force law** states that a charged particle moving in an electromagnetic
field experiences a force that depends on the values of the electric field and the
magnetic field.

**Notation:**

#. :math:`\\left[ \\vec a, \\vec b \\right]` (:code:`cross(a, b)`) is the cross product between
   :math:`\\vec a` and :math:`\\vec b`.

**Notes:**

#. This law is valid even in the relativistic case.
#. This law works only in principle because a real particle would generate its own
   electromagnetic field that would interact with the external one which would alter
   the electromagnetic force it experiences.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Lorentz_force>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

lorentz_force = clone_as_vector_symbol(symbols.force)
"""
Vector of the Lorentz :symbols:`force` exerted on the charged particle.
"""

charge = symbols.charge
"""
Value of the electric :symbols:`charge` of the test particle.
"""

electric_field = clone_as_vector_symbol(symbols.electric_field_strength)
"""
Vector of the electric field. See :symbols:`electric_field_strength`.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the particle's velocity. See :symbols:`speed`.
"""

magnetic_flux_density = clone_as_vector_symbol(symbols.magnetic_flux_density)
"""
Vector of the :symbols:`magnetic_flux_density`.
"""

law = Eq(lorentz_force, charge * (electric_field + VectorCross(velocity, magnetic_flux_density)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    charge_=charge,
    electric_field_=electric_field,
    magnetic_flux_density_=magnetic_flux_density,
    velocity_=velocity,
)
@validate_output(lorentz_force)
def calculate_lorentz_force(
    charge_: Quantity,
    electric_field_: QuantityCoordinateVector,
    magnetic_flux_density_: QuantityCoordinateVector,
    velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        charge: charge_,
        electric_field: electric_field_,
        magnetic_flux_density: magnetic_flux_density_,
        velocity: velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)
