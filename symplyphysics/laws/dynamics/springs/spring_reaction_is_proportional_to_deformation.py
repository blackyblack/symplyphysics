"""
Spring reaction is proportional to deformation
==============================================

Also called Hooke's law, it is an empirical law which states that the force needed to
extend or compress a spring by some distance is proportional with respect to the
deformation of the spring.

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

from sympy import Eq, sympify
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    Vector,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics.springs.vector import (
    spring_reaction_is_proportional_to_deformation as hookes_vector_law,)

# Also see its [vector counterpart](../vector/spring_reaction_from_deformation.py)

spring_reaction = symbols.force
"""
Restoring :symbols:`force`, or spring reaction, exerted by the spring due to the
deformation.
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

_deformation_vector = Vector([deformation])
_spring_reaction_vector_derived = hookes_vector_law.force_law(_deformation_vector)
assert len(_spring_reaction_vector_derived.components) == 1
_spring_reaction_derived = sympify(_spring_reaction_vector_derived.components[0]).subs(
    hookes_vector_law.stiffness, stiffness)
assert expr_equals(_spring_reaction_derived, law.rhs)


@validate_input(stiffness_=stiffness, deformation_=deformation)
@validate_output(spring_reaction)
def calculate_spring_reaction(stiffness_: Quantity, deformation_: Quantity) -> Quantity:
    result = law.rhs.subs({stiffness: stiffness_, deformation: deformation_})
    return Quantity(result)
