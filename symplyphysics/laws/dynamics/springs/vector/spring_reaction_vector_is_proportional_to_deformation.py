"""
Spring reaction is proportional to deformation (vector)
=======================================================

The empirical law known as **Hooke's law** states that the force needed to extend or compress
a spring by some distance scales linearly with respect to that distance. Also see the :ref:`scalar
counterpart <Spring reaction is proportional to deformation>` of this law.

**Conditions:**

#. The deformations are elastic (reversible).

#. The stiffness coefficient is a scalar.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Hooke's_law#Linear_springs>`__.
"""

from sympy import Eq
from symplyphysics import units, Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

force = clone_as_vector_symbol(symbols.force)
"""
Vector of the :symbols:`force` exerted onto the end of the spring.
"""

stiffness = symbols.stiffness
"""
:symbols:`stiffness` coefficient of the spring.
"""

deformation = clone_as_vector_symbol(symbols.distance)
"""
Vector of the deformation of the string, which is the difference between the position vector of
the end of the spring after and before the deformation. See :symbols:`distance`.
"""

law = Eq(force, -1 * stiffness * deformation)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(coefficient_=stiffness, deformation_=units.length)
@validate_output(units.force)
def calculate_force(
    coefficient_: Quantity,
    deformation_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = solve_for_vector(law, force).subs({
        stiffness: coefficient_,
        deformation: deformation_,
    })

    return QuantityCoordinateVector.from_expr(result)


@validate_input(coefficient_=stiffness, force_=units.force)
@validate_output(units.length)
def calculate_deformation(
    coefficient_: Quantity,
    force_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = solve_for_vector(law, deformation).subs({
        stiffness: coefficient_,
        force: force_,
    })

    return QuantityCoordinateVector.from_expr(result)
