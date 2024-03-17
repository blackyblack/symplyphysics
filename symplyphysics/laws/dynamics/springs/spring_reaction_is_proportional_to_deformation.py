from sympy import Eq, sympify, solve
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    print_expression,
    Vector,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics.springs.vector import (
    spring_reaction_is_proportional_to_deformation as hookes_vector_law,
)

# Description
## Also called Hooke's law, it is an empirical law which states that the force needed to extend
## or compress a spring by some distance is proportional with respect to the deformation of the
## spring.

# Law: F = -k * x
## F - restoring force (spring reaction) exerted by the spring due to the deformation
## k - proportionality coefficient, characteristic of the spring, also known as stiffness
## x - deformation of the spring (positive or negative, see Note)

# Conditions
## - `x` is small compared to the total possible deformation of the spring
## - Only applies to elastic deformations of the body, i.e. the body reverts to its initial
##   state after the removal of force or load applied onto it.

# Note
## - The spring is aligned in the positive direction of the x-axis, thus the deformation can
##   be positive, in which case the spring is stretched, or negative, in which case the spring
##   is compressed. The sign of the force indicates its direction along the x-axis.

# Also see its [vector counterpart](../vector/spring_reaction_from_deformation.py)

spring_reaction = Symbol("spring_reaction", units.force)
stiffness = Symbol("stiffness", units.force / units.length)
deformation = Symbol("deformation", units.length)

law = Eq(spring_reaction, -1 * stiffness * deformation)

# Derive current law from its vector counterpart

_deformation_vector = Vector([deformation])
_spring_reaction_vector_derived = hookes_vector_law.force_law(_deformation_vector)
assert len(_spring_reaction_vector_derived.components) == 1
_spring_reaction_derived = sympify(_spring_reaction_vector_derived.components[0]).subs(
    hookes_vector_law.stiffness, stiffness
)
assert expr_equals(_spring_reaction_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(stiffness_=stiffness, deformation_=deformation)
@validate_output(spring_reaction)
def calculate_spring_reaction(stiffness_: Quantity, deformation_: Quantity) -> Quantity:
    result = law.rhs.subs({stiffness: stiffness_, deformation: deformation_})
    return Quantity(result)
