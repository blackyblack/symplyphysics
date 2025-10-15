"""
Spring reaction is proportional to deformation
==============================================

Also called Hooke's law, it is an empirical law which states that the force needed to
extend or compress a spring by some distance is proportional with respect to the
deformation of the spring. Also see the :ref:`vector counterpart <Spring reaction is proportional
to deformation (vector)>` of this law.

**Notes:**

#. The spring is aligned in the positive direction of the :math:`x`-axis, thus the
   deformation can be positive, in which case the spring is stretched, or negative, in
   which case the spring is compressed. The sign of the force indicates its direction
   along the :math:`x`-axis.

**Conditions:**

#. :math:`x` is small compared to the total possible deformation of the spring
#. Only applies to elastic deformations of the body, i.e. the body reverts to its
   initial state after the removal of force or load applied onto it.

**Links:**

#. `Wikipedia, last formula in paragraph <https://en.wikipedia.org/wiki/Hooke%27s_law#Linear_springs>`__.
"""

from sympy import Eq
from symplyphysics import symbols, Quantity, validate_input, validate_output
from symplyphysics.laws.dynamics.springs.vector import (
    spring_reaction_vector_is_proportional_to_deformation as _hookes_vector_law,)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector, vector_equals

spring_reaction = symbols.force
"""
Restoring :symbols:`force`, or spring reaction, exerted by the spring due to the deformation. Note
that this is the *projection* of the reaction force on the displacement vector of the spring's end.
"""

stiffness = symbols.stiffness
"""
:symbols:`stiffness` of the spring.
"""

deformation = symbols.deformation
"""
:symbols:`deformation` of the spring.
"""

law = Eq(spring_reaction, -1 * stiffness * deformation)
"""
:laws:symbol::

:laws:latex::
"""

# Derive current law from its vector counterpart

_deformation_vector = CoordinateVector([deformation, 0, 0], CARTESIAN)

_spring_reaction_vector_derived = solve_for_vector(
    _hookes_vector_law.law,
    _hookes_vector_law.force,
).subs({
    _hookes_vector_law.stiffness: stiffness,
    _hookes_vector_law.deformation: _deformation_vector,
})
_spring_reaction_vector_derived = CoordinateVector.from_expr(_spring_reaction_vector_derived)

_expected_force_expr = CoordinateVector([law.rhs, 0, 0], CARTESIAN)

assert vector_equals(_spring_reaction_vector_derived, _expected_force_expr)


@validate_input(stiffness_=stiffness, deformation_=deformation)
@validate_output(spring_reaction)
def calculate_spring_reaction(stiffness_: Quantity, deformation_: Quantity) -> Quantity:
    result = law.rhs.subs({stiffness: stiffness_, deformation: deformation_})
    return Quantity(result)


# UNIQUE_LAW_ID: 250
