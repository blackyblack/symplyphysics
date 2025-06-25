"""
Force acting on dipole in non-uniform electric field
====================================================

If an electric dipole is positioned in a spatially non-uniform electric field, the forces acting
on the point charges that compose the dipole no longer cancel each other and the dipole
experiences an overall non-zero acceleration.

**Notes:**

#. A more general representation of this law, which does not require choosing an axis aligned with
   the dipole, assuming Cartesian coordinates:

   .. math::

       \\vec F = \\left( \\vec p, \\nabla \\right) \\vec E 
               = p_x \\frac{\\partial \\vec E}{\\partial x}
               + p_y \\frac{\\partial \\vec E}{\\partial y}
               + p_z \\frac{\\partial \\vec E}{\\partial z}

**Links:**

#. `Physics Bootcamp <https://www.physicsbootcamp.org/Torque-on-an-Electric-Dipole.html>`__.
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, Quantity, symbols

from symplyphysics.core.experimental.vectors import (clone_as_vector_symbol,
    clone_as_vector_function, VectorDerivative)
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

force = clone_as_vector_symbol(symbols.force)
"""
:symbols:`force` acting on the dipole.
"""

electric_dipole_moment = symbols.electric_dipole_moment
"""
Magnitude of the :symbols:`electric_dipole_moment` vector.
"""

position = symbols.position
"""
:symbols:`position` along the axis whose direction is aligned with that of the
:attr:`~electric_dipole_moment` vector.
"""

electric_field = clone_as_vector_function(symbols.electric_field_strength, (position,))
"""
Vector of the electric field as a function of :attr:`~position`. See
:symbols:`electric_field_strength`.
"""

# NOTE: we can probably make `force` a function of `position`, too, after vector integrals are
# implemented
# NOTE: `electric_dipole_moment` is independent of `position`.

law = Eq(
    force,
    electric_dipole_moment * VectorDerivative(electric_field(position), position),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    electric_dipole_moment_=electric_dipole_moment,
    position_change_=position,
    electric_field_change_=electric_field,
)
@validate_output(force)
def calculate_force(
    electric_dipole_moment_: Quantity,
    position_change_: Quantity,
    electric_field_change_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    electric_field_ = position / position_change_ * electric_field_change_

    result = law.rhs.subs({
        electric_field(position): electric_field_,
        electric_dipole_moment: electric_dipole_moment_,
    }).doit()

    return QuantityCoordinateVector.from_expr(result)
