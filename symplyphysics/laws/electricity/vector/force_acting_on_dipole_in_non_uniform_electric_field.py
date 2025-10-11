"""
Force acting on dipole in non-uniform electric field
====================================================

If an electric dipole is positioned in a spatially non-uniform electric field, the forces acting
on the point charges that compose the dipole no longer cancel each other and the dipole
experiences an overall non-zero acceleration.

**Conditions:**

#. The vector of electric dipole moment does not depent on coordinates.

**Notes:**

#.  A more general representation of this law, which does not require choosing an axis aligned with
    the dipole, assuming Cartesian coordinates:

    .. math::

        \\vec F = \\left( \\vec p, \\nabla \\right) \\vec E 
                = p_x \\frac{\\partial \\vec E}{\\partial x}
                + p_y \\frac{\\partial \\vec E}{\\partial y}
                + p_z \\frac{\\partial \\vec E}{\\partial z}

    In general orthogonal curvilinear coordinates, this equates the contraction of the *covariant
    derivative* tensor of the electric field with the electric dipole moment vector:

    .. math::

        F^i = \\left( \\nabla \\vec E : \\vec p \\right)^i = \\sum_j \\nabla_j E^i p_j,

    .. math::
    
        \\nabla_j E^i = \\partial_j E^i + \\sum_k \\Gamma^i{}_{jk} E^k.
        
    Here, superscript denotes vector components, :math:`\\partial_j` the partial derivative
    with respect to :math:`q^j`, and :math:`\\Gamma^i{}_{jk}` the **Christoffel symbols of
    the second kind** (`Wikipedia
    <https://en.wikipedia.org/wiki/Curvilinear_coordinates#Christoffel_symbols>`__). If vectors
    are represented in an *orthonormal* basis, then :math:`\\hat{\\Gamma}^i{}_{jk}` is
    used instead, related by the equation

    .. math::

        h_j h_k \\hat{\\Gamma}^i{}_{jk} = h_i \\Gamma^i{}_{jk} + \\delta_{ik} \\partial_j h_i
                                        - \\delta_{ij} \\partial_k h_j

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
    # NOTE: assuming the electric field is a linear function of position
    electric_field_ = position / position_change_ * electric_field_change_

    result = law.rhs.subs({
        electric_field(position): electric_field_,
        electric_dipole_moment: electric_dipole_moment_,
    }).doit()

    return QuantityCoordinateVector.from_expr(result)
